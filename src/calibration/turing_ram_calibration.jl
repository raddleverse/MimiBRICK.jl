using Turing
import DynamicPPL
using Distributions
using DataFrames
using CSVFiles
using LinearAlgebra
using Random
using Serialization
using Dates
using Statistics

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Runs a per-block Robust Adaptive Metropolis (RAM; Vihola 2012) calibration of the
# `sneasybrick_model` Turing model defined in `create_turing_model_sneasybrick.jl`.
#
# A single joint proposal over all ~72 free parameters is essentially never accepted, and using a
# *fixed* per-block step size (whether hand-tuned or copied from the non-Turing RAM sampler's own
# converged posterior covariance) doesn't work either: that's the *marginal* covariance
# (integrating over every other parameter), whereas a Gibbs block needs a step suited to its
# *conditional* distribution (given the other blocks' current values), which can be substantially
# narrower. There's no way to know the right fixed covariance for a Gibbs-conditional block ahead
# of time without already having a well-mixed chain.
#
# This sidesteps that entirely with an ADAPTIVE proposal per block: each block carries its own
# running Cholesky factor of its proposal covariance, which is updated after every single-block
# step using the same rule as RobustAdaptiveMetropolisSampler.jl (already a dependency of this
# package, used by the non-Turing `calibration.jl` pipeline): nudge the covariance toward whatever
# makes the *observed* per-block acceptance rate converge to `opt_α` (0.234, the theoretically
# optimal rate for high-dimensional random-walk Metropolis). Because this learns from actual
# accept/reject outcomes of that specific block, conditional on the other blocks' current state, it
# automatically finds the right conditional scale -- no need to know it in advance.
#
# This bypasses Turing's `Gibbs`/`MH` machinery entirely (custom sampling loop below) since
# `RobustAdaptiveMetropolisSampler.RAM_sample` is a standalone full-run driver, not a
# per-step/pluggable API, so its 6-line R1/R2/R3 update rule is reimplemented directly here. The
# model's log-density is evaluated via `DynamicPPL.VarInfo(model, InitFromParams(...))`, keyed by
# parameter *name* rather than position -- this also sidesteps a real pitfall from the covariance-
# matrix `MH(Σ)` approach, where the covariance is silently applied in the model's internal `~`
# declaration order rather than however you listed the block's variables.
#
# `R` (the cross-stream correlation matrix) is left on a simple LKJ(7,η) prior-redraw Metropolis
# step (not RAM-adapted): a well-scaled *local* random-walk step on a correlation matrix needs a
# proposal that respects its manifold (unit diagonal, positive-definite), which is a separate
# problem from what RAM solves here.
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

"""
    BlockState

Mutable per-block RAM state: the block's variable names (as they appear in the parameter `Dict`),
its flattened dimension, and its current Cholesky factor of the proposal covariance (adapted every
step).
"""
mutable struct BlockState
    names::Tuple
    dim::Int
    L::LowerTriangular{Float64, Matrix{Float64}} # current Cholesky factor of the proposal covariance
end

# --- Flatten/unflatten a block's variables to/from a plain vector, given the full parameter Dict. ---
function flatten_block(state::Dict{Symbol, Any}, names::Tuple)
    return vcat([state[n] isa AbstractVector ? state[n] : [state[n]] for n in names]...)
end

function unflatten_block!(state::Dict{Symbol, Any}, names::Tuple, x::Vector{Float64})
    i = 1
    for n in names
        if state[n] isa AbstractVector
            len = length(state[n])
            state[n] = x[i:(i+len-1)]
            i += len
        else
            state[n] = x[i]
            i += 1
        end
    end
    return state
end

"""
    ram_block_step!(logdensity, current_lp, state, block; opt_α=0.234, γ=2/3, iter)

Perform one Robust Adaptive Metropolis step (Vihola 2012) for a single block: propose via the
block's current Cholesky factor, accept/reject against `logdensity` (holding every other parameter
fixed at its current value in `state`), then adapt the Cholesky factor toward `opt_α`. Mutates
`state` (if accepted) and `block.L` (always). Returns `(accepted::Bool, log_density)`.
"""
function ram_block_step!(logdensity, current_lp::Float64, state::Dict{Symbol, Any}, block::BlockState;
                          opt_α::Float64 = 0.234, γ::Float64 = 2 / 3, iter::Int)
    x = flatten_block(state, block.names)
    u = randn(block.dim)
    y = x .+ block.L * u

    proposed_state = copy(state)
    unflatten_block!(proposed_state, block.names, y)
    proposed_lp = logdensity(proposed_state)

    α = min(1.0, exp(proposed_lp - current_lp))
    accepted = rand() < α
    if accepted
        for n in block.names
            state[n] = proposed_state[n]
        end
    end

    η = min(1.0, block.dim * iter^(-γ))
    M = block.L * (I + η * (α - opt_α) * (u * u') / norm(u)^2) * block.L'
    block.L = cholesky(Symmetric(M)).L

    return accepted, (accepted ? proposed_lp : current_lp)
end

# Names, in the same order used throughout this file, of the 36 scalar parameters followed by the
# 15-element `antarctic_params` vector -- matching the row/column order of both
# `calibration_initial_values_sneasybrick_04Jun2026.csv` and
# `initial_proposal_covariance_matrix_sneasybrick_04Jun2026.csv` (both package data files, both
# outputs of the same non-Turing calibration pipeline).
const SNEASYBRICK_SCALAR_PARAM_NAMES = (
    :σ_temperature, :σ_ocean_heat, :σ_glaciers, :σ_greenland, :σ_antarctic, :σ_gmsl, :σ²_white_noise_CO₂,
    :ρ_temperature, :ρ_ocean_heat, :ρ_glaciers, :ρ_greenland, :ρ_antarctic, :ρ_gmsl, :α₀_CO₂,
    :CO₂_0, :N₂O_0, :temperature_0, :ocean_heat_0, :thermal_s₀, :greenland_v₀, :glaciers_v₀, :glaciers_s₀, :antarctic_s₀,
    :Q10, :CO₂_fertilization, :CO₂_diffusivity,
    :heat_diffusivity, :rf_scale_aerosol, :ECS,
    :thermal_α,
    :greenland_a, :greenland_b, :greenland_α, :greenland_β,
    :glaciers_β₀, :glaciers_n,
)

const SNEASYBRICK_GIBBS_BLOCKS = (
    sneasy    = (:CO₂_0, :N₂O_0, :temperature_0, :ocean_heat_0, :Q10, :CO₂_fertilization, :CO₂_diffusivity, :heat_diffusivity, :rf_scale_aerosol, :ECS),
    glaciers  = (:glaciers_v₀, :glaciers_s₀, :glaciers_β₀, :glaciers_n),
    thermal   = (:thermal_s₀, :thermal_α),
    antarctic = (:antarctic_s₀, :antarctic_params),
    greenland = (:greenland_v₀, :greenland_a, :greenland_b, :greenland_α, :greenland_β),
    autocov   = (:σ_temperature, :σ_ocean_heat, :σ_glaciers, :σ_greenland, :σ_antarctic, :σ_gmsl,
                 :σ²_white_noise_CO₂, :ρ_temperature, :ρ_ocean_heat, :ρ_glaciers, :ρ_greenland,
                 :ρ_antarctic, :ρ_gmsl, :α₀_CO₂),
)

# Maps a Gibbs block's variable names to their integer position(s) in
# `SNEASYBRICK_SCALAR_PARAM_NAMES ∪ antarctic_params[1:15]` (positions 1:36 and 37:51
# respectively), so per-block warm-start covariances can be sliced directly out of the bundled
# `initial_proposal_covariance_matrix_sneasybrick_04Jun2026.csv`.
function block_param_indices(names::Tuple)
    param_index = Dict(name => i for (i, name) in enumerate(SNEASYBRICK_SCALAR_PARAM_NAMES))
    idx = Int[]
    for n in names
        if n === :antarctic_params
            append!(idx, 37:51)
        else
            push!(idx, param_index[n])
        end
    end
    return idx
end

# The same 51 quantities as `SNEASYBRICK_SCALAR_PARAM_NAMES ∪ antarctic_params[1:15]`, in the same
# order, but named as they appear in the package's established calibration-output convention (e.g.
# `parameters_subsample_sneasybrick.csv`, consumed by `run_projections`). Used to translate this
# sampler's internal parameter Dict into that convention.
const SNEASYBRICK_EXTERNAL_PARAM_NAMES = (
    "sd_temp","sd_ocean_heat","sd_glaciers","sd_greenland","sd_antarctic","sd_gmsl",
    "sigma_whitenoise_co2","rho_temperature","rho_ocean_heat","rho_glaciers","rho_greenland",
    "rho_antarctic","rho_gmsl","alpha0_CO2","CO2_0","N2O_0","temperature_0","ocean_heat_0",
    "thermal_s0","greenland_v0","glaciers_v0","glaciers_s0","antarctic_s0","Q10",
    "CO2_fertilization","CO2_diffusivity","heat_diffusivity","rf_scale_aerosol","climate_sensitivity",
    "thermal_alpha","greenland_a","greenland_b","greenland_alpha","greenland_beta","glaciers_beta0",
    "glaciers_n","anto_alpha","anto_beta","antarctic_gamma","antarctic_alpha","antarctic_mu",
    "antarctic_nu","antarctic_precip0","antarctic_kappa","antarctic_flow0","antarctic_runoff_height0",
    "antarctic_c","antarctic_bed_height0","antarctic_slope","antarctic_lambda","antarctic_temp_threshold",
)

"""
    run_per_block_ram_calibration(n_samples; kwargs...)

Sample the SNEASY-BRICK Turing model (`sneasybrick_model`) with a per-block Robust Adaptive
Metropolis sampler (see module docstring above): SNEASY, glaciers, thermal expansion, Antarctic ice
sheet, Greenland ice sheet, and the statistical auto-covariance parameters each get their own
adaptively-tuned multivariate proposal; `R` (cross-covariance) is a simple LKJ(7,η) prior-redraw
Metropolis step.

Function Arguments:

    n_samples                = Number of outer sweeps (one RAM step per block, per sweep) to draw.
    output_dir                = Directory the chain (and a timing summary) are saved to.
    calibration_start_year   = First year to run the model (default 1850).
    calibration_end_year     = Final year of the calibration (default 2017).
    joint_antarctic_prior    = TRUE/FALSE: joint normal (TRUE) vs. marginal kernel density (FALSE, default)
                                Antarctic ice sheet prior.
    uniform_ECS               = TRUE/FALSE: uniform (TRUE) vs. paleo-informed Cauchy (FALSE, default) ECS prior.
    eta                       = Fixed `LKJ(7, η)` concentration for the `R` prior. `nothing` (default) chooses `η`
                                from the empirically estimated cross-stream correlation.
    reference_posterior_path  = Optional path to a CSV of posterior parameter draws (e.g. a previous
                                `run_calibration` subsample, such as `parameters_subsample_sneasybrick.csv`)
                                used to build each block's initial proposal covariance. Defaults to `nothing`,
                                which instead seeds every block from the package's own bundled
                                `initial_proposal_covariance_matrix_sneasybrick_04Jun2026.csv`. Either way this
                                is only a warm start -- RAM adapts away a poorly-scaled starting covariance
                                over the run.
    calibration_data_dir      = Directory containing the package's calibration data files. Defaults to the
                                bundled `data/calibration_data` directory.
    seed                      = RNG seed for reproducibility.

Returns `(chain_df, elapsed_seconds, acceptance_rates)`, where `chain_df` is a `DataFrame` with one
column per (scalar or vector-element) parameter and one row per sweep, and `acceptance_rates` is a
`Dict{Symbol,Float64}` of each block's overall acceptance rate.
"""
function run_per_block_ram_calibration(n_samples::Int;
                                        output_dir::String,
                                        calibration_start_year::Int = 1850,
                                        calibration_end_year::Int = 2017,
                                        joint_antarctic_prior::Bool = false,
                                        uniform_ECS::Bool = false,
                                        eta::Union{Real, Nothing} = nothing,
                                        reference_posterior_path::Union{String, Nothing} = nothing,
                                        calibration_data_dir::Union{String, Nothing} = nothing,
                                        seed::Int = 2024)

    mkpath(output_dir)

    if isnothing(calibration_data_dir)
        calibration_data_dir = joinpath(@__DIR__, "..", "..", "data", "calibration_data")
    end

    run_mymodel! = construct_run_sneasybrick(calibration_start_year, calibration_end_year)
    model, empirical_R, eta = construct_sneasybrick_turing_model(run_mymodel!;
        model_start_year = calibration_start_year, calibration_end_year = calibration_end_year,
        joint_antarctic_prior = joint_antarctic_prior, uniform_ECS = uniform_ECS, eta = eta,
        calibration_data_dir = calibration_data_dir)
    println("Using η = $eta for the R ~ LKJ(7, η) prior.")

    #---------------------------------------------------------------------------
    # Seed the chain from the package's known-good calibration starting values (raw prior draws are
    # essentially always numerically invalid for this model).
    #---------------------------------------------------------------------------
    initvals = DataFrame(load(joinpath(calibration_data_dir, "calibration_initial_values_sneasybrick_04Jun2026.csv"))).parameter_values
    state = Dict{Symbol, Any}(zip(SNEASYBRICK_SCALAR_PARAM_NAMES, initvals[1:36]))
    state[:antarctic_params] = Vector{Float64}(initvals[37:51])
    # Start `R` at the identity (no cross-stream correlation) rather than the empirical estimate:
    # `empirical_R` combined with the calibration_initial_values σ/ρ makes Σ_total non-PD at this
    # exact point, giving -Inf at sweep 1 and forcing every *other* block to reject until R happens
    # to redraw into a compatible state. Identity reliably gives a finite starting log-density, so
    # every block's RAM adaptation can start immediately instead of waiting on R.
    state[:R] = Matrix{Float64}(I, 7, 7)

    #---------------------------------------------------------------------------
    # Initial per-block proposal covariances (a warm start only -- RAM adapts away a poorly-scaled
    # guess over the run). Default: slice the package's own bundled covariance matrix, which shares
    # row/column order and units (including log-space `antarctic_precip0`) with
    # `calibration_initial_values_sneasybrick_04Jun2026.csv`. If the caller supplies
    # `reference_posterior_path` (e.g. a previous calibration's posterior subsample), use its
    # empirical covariance instead; that file stores `antarctic_precip0` in raw physical units, so
    # it's log-transformed first to match the model's internal log-space convention.
    #---------------------------------------------------------------------------
    if reference_posterior_path === nothing
        full_cov = Matrix(DataFrame(load(joinpath(calibration_data_dir, "initial_proposal_covariance_matrix_sneasybrick_04Jun2026.csv"))))
    else
        ram_df = DataFrame(load(reference_posterior_path))
        ram_df.antarctic_precip0 = log.(ram_df.antarctic_precip0)
        full_cov = Matrix(cov(Matrix(Float64.(ram_df[!, collect(SNEASYBRICK_EXTERNAL_PARAM_NAMES)]))))
    end
    initial_cov(idx) = Symmetric((2.38^2 / length(idx)) .* full_cov[idx, idx])

    blocks = NamedTuple(
        block_name => BlockState(names, length(flatten_block(state, names)), cholesky(initial_cov(block_param_indices(names))).L)
        for (block_name, names) in pairs(SNEASYBRICK_GIBBS_BLOCKS)
    )

    #---------------------------------------------------------------------------
    # Log-density evaluator: rebuilds a VarInfo from the full named-parameter Dict every call (so
    # this pays the same per-evaluation cost -- one SNEASY-BRICK simulation -- as a Gibbs/MH
    # approach; no efficiency lost by not using Turing's own Gibbs machinery here).
    #---------------------------------------------------------------------------
    function logdensity(s::Dict{Symbol, Any})
        nt = NamedTuple(s)
        vi = DynamicPPL.VarInfo(model, DynamicPPL.InitFromParams(nt, nothing))
        lp = DynamicPPL.getlogp(vi)
        return lp.logprior + lp.logjac + lp.loglikelihood
    end

    #---------------------------------------------------------------------------
    # R: a simple (non-adaptive) LKJ(7,η) prior-redraw Metropolis step.
    #---------------------------------------------------------------------------
    r_prior = LKJ(7, eta)
    function r_step!(current_lp::Float64, state::Dict{Symbol, Any})
        proposed_state = copy(state)
        R_proposed = rand(r_prior)
        proposed_state[:R] = R_proposed
        proposed_lp = logdensity(proposed_state)
        α = min(1.0, exp(proposed_lp - current_lp))
        if rand() < α
            state[:R] = R_proposed
            return true, proposed_lp
        end
        return false, current_lp
    end

    #---------------------------------------------------------------------------
    # Sampling loop: one RAM step per block, plus one R redraw, per outer sweep.
    #---------------------------------------------------------------------------
    Random.seed!(seed)
    current_lp = logdensity(state)
    n_accept = Dict{Symbol, Int}(name => 0 for name in keys(blocks))
    n_accept[:R] = 0

    # Storage: one column per scalar parameter, plus antarctic_params[1..15] and R_i_j, plus logjoint.
    scalar_cols = Dict{Symbol, Vector{Float64}}(name => Vector{Float64}(undef, n_samples) for name in SNEASYBRICK_SCALAR_PARAM_NAMES)
    antarctic_cols = [Vector{Float64}(undef, n_samples) for _ in 1:15]
    r_cols = [Vector{Float64}(undef, n_samples) for _ in 1:49]
    logjoint_col = Vector{Float64}(undef, n_samples)

    println("Sampling $n_samples sweeps with per-block RAM (SNEASY / glaciers / thermal / Antarctic / Greenland / auto-covariance / cross-covariance)...")
    elapsed = @elapsed for iter in 1:n_samples
        for (block_name, block) in pairs(blocks)
            accepted, current_lp = ram_block_step!(logdensity, current_lp, state, block; iter = iter)
            n_accept[block_name] += accepted
        end
        accepted_r, current_lp = r_step!(current_lp, state)
        n_accept[:R] += accepted_r

        for name in SNEASYBRICK_SCALAR_PARAM_NAMES
            scalar_cols[name][iter] = state[name]
        end
        for i in 1:15
            antarctic_cols[i][iter] = state[:antarctic_params][i]
        end
        for i in 1:7, j in 1:7
            r_cols[(i-1)*7 + j][iter] = state[:R][i, j]
        end
        logjoint_col[iter] = current_lp

        if iter % max(1, n_samples ÷ 20) == 0
            println("  sweep $iter/$n_samples  (", round(100*iter/n_samples, digits=0), "%)  logjoint=$(round(current_lp,digits=1))")
        end
    end
    println("Done: $n_samples sweeps in $(round(elapsed, digits=1))s ($(round(elapsed/n_samples, digits=3)) s/sweep)")

    acceptance_rates = Dict(name => n_accept[name] / n_samples for name in keys(n_accept))
    println("Acceptance rates: ", acceptance_rates)

    #---------------------------------------------------------------------------
    # Assemble and save the chain as a DataFrame.
    #---------------------------------------------------------------------------
    df = DataFrame()
    for name in SNEASYBRICK_SCALAR_PARAM_NAMES
        df[!, string(name)] = scalar_cols[name]
    end
    for i in 1:15
        df[!, "antarctic_params[$i]"] = antarctic_cols[i]
    end
    for i in 1:7, j in 1:7
        df[!, "R_$(i)_$(j)"] = r_cols[(i-1)*7 + j]
    end
    df[!, "logjoint"] = logjoint_col

    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    basename_ = "sneasybrick_ram_chain_n$(n_samples)_$(timestamp)"

    csv_path = joinpath(output_dir, basename_ * ".csv")
    save(csv_path, df)
    println("Saved chain as CSV to $csv_path")

    chain_path = joinpath(output_dir, basename_ * ".jls")
    serialize(chain_path, df)
    println("Saved chain to $chain_path")

    timing_path = joinpath(output_dir, basename_ * "_timing.csv")
    save(timing_path, DataFrame(n_samples = [n_samples], elapsed_seconds = [elapsed], seconds_per_sweep = [elapsed / n_samples]))
    println("Saved timing summary to $timing_path")

    return df, elapsed, acceptance_rates
end
