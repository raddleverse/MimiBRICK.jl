using MimiBRICK
using Turing
using DynamicPPL
using AdvancedMH
using Distributions
using DataFrames
using CSVFiles
using LinearAlgebra
using Random
using Serialization
using Dates

##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
## Runs and benchmarks a Turing calibration of SNEASY-BRICK (see
## src/calibration/create_log_posteriors/create_turing_model_sneasybrick.jl), saving the
## resulting MCMC chain to disk.
##
## Uses Turing's gradient-free Metropolis-Hastings sampler (`MH()`), not NUTS: NUTS's default
## (ForwardDiff) autodiff backend cannot be used on this model because SNEASY-BRICK is run
## through a `Mimi.ModelInstance` whose scalar parameters are stored as concretely-typed
## `Mimi.ScalarModelParameter{Float64}` (fixed at `Mimi.build()` time), so `ForwardDiff.Dual`
## numbers can't be pushed through it. See examples/turing_nuts_calibration_test.jl for a
## reproduction of that failure. This is a Mimi.jl architecture limitation, not something
## fixable in this repo without a much heavier per-evaluation cost (rebuilding the whole Mimi
## model on every gradient call).
##
## Also NOT using bare `MH()`: with no arguments, Turing's `MH()` proposes every parameter
## fresh from its *prior* each iteration (an independence sampler), rather than taking local
## steps from the current state. Since raw prior draws are (as above) essentially always
## numerically invalid for this model, that gives a ~0% acceptance rate and a chain that never
## leaves its starting point.
##
## And NOT a single joint random-walk proposal over every parameter either: Turing's per-symbol
## `MH(:a => proposal_a, :b => proposal_b, ...)` treats ALL named parameters as ONE combined
## proposal, accepted/rejected as a whole. Across this model's free dimensions, even
## individually well-scaled per-parameter steps combine into an almost-never-accepted joint
## move (empirically ~0-1% acceptance). Instead this uses `Gibbs` to split the parameters into
## smaller blocks, each updated (proposed + accepted/rejected) in its own step, with hand-tuned
## random-walk step sizes (~2% of each parameter's prior range/scale; the naive 5-10% guess that
## works fine for a *single* parameter is still far too large once several parameters have to move
## together and jointly land in a good region):
##   - SNEASY (climate/carbon-cycle)
##   - BRICK, split by sea level component rather than one single block, since each component's
##     parameters only affect that component's own modeled contribution (glaciers, thermal
##     expansion, Antarctica, Greenland) -- splitting them lets each move independently rather
##     than needing every component to jointly land in a good region at once:
##       - Glaciers & small ice caps
##       - Thermal expansion
##       - Antarctic ice sheet (including the 15 paleo-informed ice sheet parameters)
##       - Greenland ice sheet
##   - the statistical auto-covariance parameters (the AR(1)/CAR(1) σ's and ρ's/α₀)
##
## The cross-stream correlation matrix `R` is sampled too (`R ~ LKJ(7, η)`), in its own Gibbs
## block, on Turing's default (prior/LKJ-redraw) proposal -- a well-scaled *local* proposal on a
## correlation matrix needs one that respects its manifold (unit diagonal, positive-definite),
## which AdvancedMH's generic `RandomWalkProposal` does not provide out of the box. `η` itself is
## fixed (not sampled) and chosen -- see create_turing_model_sneasybrick.jl's
## `empirical_cross_correlation` and `choose_lkj_eta` -- to be informed by the empirically
## observed cross-stream correlation strength, rather than defaulting to an uninformative `η=1`.
## `R`'s block is initialized at the empirical estimate (a much better starting point than the
## identity matrix).
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------

"""
    run_mh_calibration(n_samples::Int; kwargs...)

Sample the SNEASY-BRICK Turing model with a 7-block Gibbs sampler (SNEASY / glaciers / thermal
expansion / Antarctic ice sheet / Greenland ice sheet / auto-covariance / cross-covariance, each
a within-block `MH()`), benchmark the run, and save the resulting chain to `output_dir`.

Function Arguments:

    n_samples                = Number of MH samples to draw.
    calibration_start_year   = First year to run the model (default 1850).
    calibration_end_year     = Final year of the calibration (default 2017).
    joint_antarctic_prior    = TRUE/FALSE: joint normal (TRUE) vs. marginal kernel density (FALSE, default)
                                Antarctic ice sheet prior.
    uniform_ECS               = TRUE/FALSE: uniform (TRUE) vs. paleo-informed Cauchy (FALSE, default) ECS prior.
    eta                       = Fixed `LKJ(7, η)` concentration for the `R` prior. `nothing` (default) chooses `η`
                                from the empirically estimated cross-stream correlation; pass a specific value
                                (e.g. for a sensitivity analysis -- see turing_mh_R_sensitivity_analysis.jl) to
                                override that choice.
    output_dir                = Directory the chain (and a timing summary) are saved to.
    seed                      = RNG seed for reproducibility.

Returns `(chain, elapsed_seconds, eta)`.
"""
function run_mh_calibration(n_samples::Int;
                             calibration_start_year::Int = 1850,
                             calibration_end_year::Int = 2017,
                             joint_antarctic_prior::Bool = false,
                             uniform_ECS::Bool = false,
                             eta::Union{Real, Nothing} = nothing,
                             output_dir::String = joinpath(@__DIR__, "results"),
                             seed::Int = 2024)

    mkpath(output_dir)

    run_mymodel! = MimiBRICK.construct_run_sneasybrick(calibration_start_year, calibration_end_year)
    model, empirical_R, eta = MimiBRICK.construct_sneasybrick_turing_model(run_mymodel!;
        model_start_year = calibration_start_year, calibration_end_year = calibration_end_year,
        joint_antarctic_prior = joint_antarctic_prior, uniform_ECS = uniform_ECS, eta = eta)
    println("Using η = $eta for the R ~ LKJ(7, η) prior.")

    #---------------------------------------------------------------------------
    # Seed the chain from the package's known-good calibration starting values. Raw draws
    # from the priors are essentially always numerically invalid for this model (this is a
    # pre-existing property of SNEASY-BRICK, not specific to MH or Turing -- `run_calibration`'s
    # own `start_from_priors=true` option is explicitly unimplemented for the same reason).
    #---------------------------------------------------------------------------
    calibration_data_dir = joinpath(pkgdir(MimiBRICK), "data", "calibration_data")
    initvals = DataFrame(load(joinpath(calibration_data_dir, "calibration_initial_values_sneasybrick_04Jun2026.csv"))).parameter_values

    # Same order as `p` in `sneasybrick_model` / the initial-values CSV (36 scalar parameters,
    # then the 15 Antarctic ice sheet parameters as a block).
    scalar_param_names = (
        :σ_temperature, :σ_ocean_heat, :σ_glaciers, :σ_greenland, :σ_antarctic, :σ_gmsl, :σ²_white_noise_CO₂,
        :ρ_temperature, :ρ_ocean_heat, :ρ_glaciers, :ρ_greenland, :ρ_antarctic, :ρ_gmsl, :α₀_CO₂,
        :CO₂_0, :N₂O_0, :temperature_0, :ocean_heat_0, :thermal_s₀, :greenland_v₀, :glaciers_v₀, :glaciers_s₀, :antarctic_s₀,
        :Q10, :CO₂_fertilization, :CO₂_diffusivity,
        :heat_diffusivity, :rf_scale_aerosol, :ECS,
        :thermal_α,
        :greenland_a, :greenland_b, :greenland_α, :greenland_β,
        :glaciers_β₀, :glaciers_n,
    )

    initial_params = merge(
        NamedTuple(zip(scalar_param_names, initvals[1:36])),
        (antarctic_params = Vector{Float64}(initvals[37:51]), R = empirical_R), # a much better start than R = identity.
    )
    init_strategy = DynamicPPL.InitFromParams(initial_params, nothing)

    #---------------------------------------------------------------------------
    # Random-walk step sizes per scalar parameter, expressed as a base value (~5-10% of each
    # parameter's prior range/scale -- see the priors declared in `sneasybrick_model`,
    # create_turing_model_sneasybrick.jl) times a *per-block* scale factor.
    #
    # A single shared scale factor (originally 0.02 for every block) doesn't work well once blocks
    # have very different dimensions: an ESS/wall-clock-time sweep per block (holding every other
    # block fixed) found that low-dimensional blocks (glaciers: 4 dims, thermal: 2, greenland: 5)
    # actually do *better* with a *larger* scale than 0.02 (they were mixing fine, but taking
    # unnecessarily small, slow-to-explore steps), while higher-dimensional blocks (SNEASY: 10,
    # auto-covariance: 14, Antarctic: 16) need a *smaller* scale to keep the jointly-proposed move
    # likely enough to be accepted at all. These single-run ESS estimates are noisy (300-500
    # samples each), so treat them as directional, not exact optima.
    #---------------------------------------------------------------------------
    block_step_scales = Dict(
        :sneasy    => 0.002,  # 10 dims
        :glaciers  => 0.2,    #  4 dims
        :thermal   => 0.4,    #  2 dims
        :antarctic => 0.005,  # 16 dims (1 + the 15 paleo-informed ice sheet parameters)
        :greenland => 0.05,   #  5 dims
        :autocov   => 0.001,  # 14 dims
    )
    base_step_sizes = Dict(
        :σ_temperature => 0.02, :σ_ocean_heat => 0.4, :σ_glaciers => 0.00015, :σ_greenland => 0.0002,
        :σ_antarctic => 0.0063, :σ_gmsl => 0.005, :σ²_white_noise_CO₂ => 20.0,
        :ρ_temperature => 0.2, :ρ_ocean_heat => 0.2, :ρ_glaciers => 0.2, :ρ_greenland => 0.2,
        :ρ_antarctic => 0.2, :ρ_gmsl => 0.05, :α₀_CO₂ => 1.15,
        :CO₂_0 => 0.6, :N₂O_0 => 1.8, :temperature_0 => 0.06, :ocean_heat_0 => 10.0, :thermal_s₀ => 0.01,
        :greenland_v₀ => 0.04, :glaciers_v₀ => 0.022, :glaciers_s₀ => 0.013, :antarctic_s₀ => 0.01,
        :Q10 => 0.4, :CO₂_fertilization => 0.1, :CO₂_diffusivity => 2.0,
        :heat_diffusivity => 0.3, :rf_scale_aerosol => 0.3, :ECS => 0.5,
        :thermal_α => 0.025,
        :greenland_a => 0.4, :greenland_b => 0.3, :greenland_α => 0.0001, :greenland_β => 0.0001,
        :glaciers_β₀ => 0.004, :glaciers_n => 0.045,
    )
    step_size(block, name) = block_step_scales[block] * base_step_sizes[name]

    #---------------------------------------------------------------------------
    # Block 1: SNEASY (climate + carbon cycle).
    #---------------------------------------------------------------------------
    sneasy_vars = (:CO₂_0, :N₂O_0, :temperature_0, :ocean_heat_0, :Q10, :CO₂_fertilization,
                   :CO₂_diffusivity, :heat_diffusivity, :rf_scale_aerosol, :ECS)
    sneasy_sampler = MH([name => AdvancedMH.RandomWalkProposal(Normal(0, step_size(:sneasy, name))) for name in sneasy_vars]...)

    #---------------------------------------------------------------------------
    # BRICK, split into one block per sea level component rather than one single combined block:
    # each component's parameters only affect that component's own modeled output, so there's no
    # need to make them all move together (and, as with SNEASY vs. BRICK above, jointly landing
    # in a good region gets harder the more parameters/components have to move at once).
    #---------------------------------------------------------------------------

    # Block 2: glaciers & small ice caps.
    glaciers_vars = (:glaciers_v₀, :glaciers_s₀, :glaciers_β₀, :glaciers_n)
    glaciers_sampler = MH([name => AdvancedMH.RandomWalkProposal(Normal(0, step_size(:glaciers, name))) for name in glaciers_vars]...)

    # Block 3: thermal expansion.
    thermal_vars = (:thermal_s₀, :thermal_α)
    thermal_sampler = MH([name => AdvancedMH.RandomWalkProposal(Normal(0, step_size(:thermal, name))) for name in thermal_vars]...)

    # Block 4: Antarctic ice sheet (including the 15 paleo-informed ice sheet parameters).
    antarctic_step = block_step_scales[:antarctic] .* (model.args.antarctic_upper_bound .- model.args.antarctic_lower_bound)
    antarctic_sampler = MH(
        :antarctic_s₀ => AdvancedMH.RandomWalkProposal(Normal(0, step_size(:antarctic, :antarctic_s₀))),
        :antarctic_params => AdvancedMH.RandomWalkProposal(MvNormal(zeros(15), Diagonal(antarctic_step .^ 2))),
    )

    # Block 5: Greenland ice sheet.
    greenland_vars = (:greenland_v₀, :greenland_a, :greenland_b, :greenland_α, :greenland_β)
    greenland_sampler = MH([name => AdvancedMH.RandomWalkProposal(Normal(0, step_size(:greenland, name))) for name in greenland_vars]...)

    #---------------------------------------------------------------------------
    # Block 6: auto-covariance (the AR(1)/CAR(1) statistical nuisance parameters).
    #---------------------------------------------------------------------------
    autocov_vars = (:σ_temperature, :σ_ocean_heat, :σ_glaciers, :σ_greenland, :σ_antarctic, :σ_gmsl,
                    :σ²_white_noise_CO₂, :ρ_temperature, :ρ_ocean_heat, :ρ_glaciers, :ρ_greenland,
                    :ρ_antarctic, :ρ_gmsl, :α₀_CO₂)
    autocov_sampler = MH([name => AdvancedMH.RandomWalkProposal(Normal(0, step_size(:autocov, name))) for name in autocov_vars]...)

    #---------------------------------------------------------------------------
    # Block 7: cross-covariance (`R`, the 7x7 cross-stream correlation matrix). Left on Turing's
    # default (prior/LKJ) proposal -- see the module docstring for why.
    #---------------------------------------------------------------------------
    gibbs_sampler = Gibbs(
        sneasy_vars => sneasy_sampler,
        glaciers_vars => glaciers_sampler,
        thermal_vars => thermal_sampler,
        (:antarctic_s₀, :antarctic_params) => antarctic_sampler,
        greenland_vars => greenland_sampler,
        autocov_vars => autocov_sampler,
        (:R,) => MH(),
    )

    #---------------------------------------------------------------------------
    # Benchmark: a short timed pre-run gives a throughput estimate (and projected total time
    # for `n_samples`) before committing to the full, possibly long, run below. The first call
    # to `sample()` in a Julia session pays a one-time JIT/precompilation cost that has nothing
    # to do with steady-state sampling speed (and can dwarf it), so that cost is paid off with a
    # small *untimed* warm-up before the *timed* benchmark batch.
    #---------------------------------------------------------------------------
    println("Warming up (one-time JIT/precompilation, untimed)...")
    Random.seed!(seed)
    sample(model, gibbs_sampler, 2; initial_params = init_strategy, progress = false)

    n_bench = min(20, n_samples)
    println("Benchmarking throughput with $n_bench MH sample(s)...")
    Random.seed!(seed)
    bench_chain = nothing
    bench_elapsed = @elapsed bench_chain = sample(model, gibbs_sampler, n_bench; initial_params = init_strategy, progress = false)
    bench_rate = bench_elapsed / n_bench
    println("  $n_bench samples in $(round(bench_elapsed, digits=1))s => $(round(bench_rate, digits=3)) s/sample")
    println("  projected time for $n_samples samples: $(round(bench_rate * n_samples / 60, digits=1)) min ($(round(bench_rate * n_samples / 3600, digits=2)) hr)")

    # Sanity check: catch a degenerate (frozen / non-exploring) chain -- in one representative
    # variable from each of the seven Gibbs blocks -- before committing to the full run below.
    sanity_check_vars = ("ECS", "glaciers_n", "thermal_α", "antarctic_params[15]", "greenland_a", "σ_temperature", "R[2, 1]")
    print_sanity_check(chn, n) = for name in sanity_check_vars
        n_unique = length(unique(vec(Array(chn[:, name, :]))))
        println("  $name: $n_unique / $n unique values", n_unique == 1 ? "  <-- WARNING: not moving!" : "")
    end
    println("  sanity check (benchmark batch):")
    print_sanity_check(bench_chain, n_bench)

    #---------------------------------------------------------------------------
    # The actual, timed run.
    #---------------------------------------------------------------------------
    println("Sampling $n_samples draws with the SNEASY / glaciers / thermal / Antarctic / Greenland / auto-covariance / cross-covariance Gibbs blocks...")
    Random.seed!(seed)
    elapsed = @elapsed chain = sample(model, gibbs_sampler, n_samples; initial_params = init_strategy)
    rate = elapsed / n_samples
    println("Done: $n_samples samples in $(round(elapsed, digits=1))s ($(round(rate, digits=3)) s/sample, $(round(1/rate, digits=3)) samples/s)")

    println("Sanity check (full chain):")
    print_sanity_check(chain, n_samples)

    #---------------------------------------------------------------------------
    # Save the chain: a full-fidelity Julia serialization (reload with `Serialization.deserialize`
    # for further Turing/MCMCChains analysis) and a portable CSV (parameter draws + log-density
    # columns), matching this package's existing convention of saving calibration output as CSV
    # (see src/calibration/calibration.jl).
    #---------------------------------------------------------------------------
    timestamp  = Dates.format(now(), "yyyymmdd_HHMMSS")
    basename_  = "sneasybrick_mh_chain_n$(n_samples)_$(timestamp)"

    chain_path = joinpath(output_dir, basename_ * ".jls")
    serialize(chain_path, chain)
    println("Saved full chain to $chain_path")

    csv_path = joinpath(output_dir, basename_ * ".csv")
    save(csv_path, DataFrame(chain))
    println("Saved chain as CSV to $csv_path")

    timing_path = joinpath(output_dir, basename_ * "_timing.csv")
    save(timing_path, DataFrame(n_samples = [n_samples], elapsed_seconds = [elapsed], seconds_per_sample = [rate], samples_per_second = [1 / rate]))
    println("Saved timing summary to $timing_path")

    return chain, elapsed, eta
end

#-------------------------------------------------------------------------------
# Default: 10,000 samples. Override either by passing a sample count on the command line
# (`julia --project=. examples/turing_mh_calibration_benchmark.jl 500`) or by calling
# `run_mh_calibration(N)` directly from the REPL with any other keyword arguments.
#-------------------------------------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    n_samples = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 10_000
    run_mh_calibration(n_samples)
end
