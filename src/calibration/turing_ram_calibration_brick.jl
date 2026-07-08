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
# Runs a per-block Robust Adaptive Metropolis (RAM; Vihola 2012) calibration of the `brick_model`
# Turing model defined in `create_turing_model_brick.jl` -- the standalone-BRICK analogue of
# `run_per_block_ram_calibration` (see that function's docstring in `turing_ram_calibration.jl` for
# the full design rationale: fixed per-block Gibbs step sizes don't work because a block's
# *conditional* distribution, given the other blocks' current values, can be much narrower than its
# *marginal* one, so each block instead carries its own adaptively-tuned proposal covariance).
#
# BRICK has no SNEASY block (it has no carbon-cycle/climate parameters -- temperature and ocean
# heat are exogenous forcing here, not calibration targets), so this reuses `BlockState`,
# `flatten_block`, `unflatten_block!`, and `ram_block_step!` from `turing_ram_calibration.jl`
# as-is (they're generic, not SNEASY-BRICK-specific) with four Gibbs blocks (glaciers, thermal
# expansion, Antarctic ice sheet, Greenland ice sheet) plus the statistical auto-covariance
# parameters for BRICK's own four correlated streams, and a simple LKJ(4,η) prior-redraw
# Metropolis step for the cross-stream correlation matrix `R`.
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Names, in the same order used throughout this file, of the 20 scalar parameters followed by the
# 15-element `antarctic_params` vector -- matching the row/column order of both
# `calibration_initial_values_brick_04Jun2026.csv` and
# `initial_proposal_covariance_matrix_brick_04Jun2026.csv` (both package data files, both outputs
# of the same non-Turing calibration pipeline).
const BRICK_SCALAR_PARAM_NAMES = (
    :σ_glaciers, :σ_greenland, :σ_antarctic, :σ_gmsl,
    :ρ_glaciers, :ρ_greenland, :ρ_antarctic, :ρ_gmsl,
    :thermal_s₀, :greenland_v₀, :glaciers_v₀, :glaciers_s₀, :antarctic_s₀, :thermal_α,
    :greenland_a, :greenland_b, :greenland_α, :greenland_β,
    :glaciers_β₀, :glaciers_n,
)

const BRICK_GIBBS_BLOCKS = (
    glaciers  = (:glaciers_v₀, :glaciers_s₀, :glaciers_β₀, :glaciers_n),
    thermal   = (:thermal_s₀, :thermal_α),
    antarctic = (:antarctic_s₀, :antarctic_params),
    greenland = (:greenland_v₀, :greenland_a, :greenland_b, :greenland_α, :greenland_β),
    autocov   = (:σ_glaciers, :σ_greenland, :σ_antarctic, :σ_gmsl, :ρ_glaciers, :ρ_greenland, :ρ_antarctic, :ρ_gmsl),
)

# Maps a Gibbs block's variable names to their integer position(s) in
# `BRICK_SCALAR_PARAM_NAMES ∪ antarctic_params[1:15]` (positions 1:20 and 21:35 respectively), so
# per-block warm-start covariances can be sliced directly out of the bundled
# `initial_proposal_covariance_matrix_brick_04Jun2026.csv`.
function block_param_indices_brick(names::Tuple)
    param_index = Dict(name => i for (i, name) in enumerate(BRICK_SCALAR_PARAM_NAMES))
    idx = Int[]
    for n in names
        if n === :antarctic_params
            append!(idx, 21:35)
        else
            push!(idx, param_index[n])
        end
    end
    return idx
end

# The same 35 quantities as `BRICK_SCALAR_PARAM_NAMES ∪ antarctic_params[1:15]`, in the same order,
# but named as they appear in the package's established calibration-output convention (e.g.
# `parameters_subsample_brick.csv`, consumed by `run_projections`). Used to translate this
# sampler's internal parameter Dict into that convention.
const BRICK_EXTERNAL_PARAM_NAMES = (
    "sd_glaciers", "sd_greenland", "sd_antarctic", "sd_gmsl",
    "rho_glaciers", "rho_greenland", "rho_antarctic", "rho_gmsl",
    "thermal_s0", "greenland_v0", "glaciers_v0", "glaciers_s0", "antarctic_s0", "thermal_alpha",
    "greenland_a", "greenland_b", "greenland_alpha", "greenland_beta",
    "glaciers_beta0", "glaciers_n",
    "anto_alpha", "anto_beta", "antarctic_gamma", "antarctic_alpha", "antarctic_mu",
    "antarctic_nu", "antarctic_precip0", "antarctic_kappa", "antarctic_flow0", "antarctic_runoff_height0",
    "antarctic_c", "antarctic_bed_height0", "antarctic_slope", "antarctic_lambda", "antarctic_temp_threshold",
)

"""
    run_per_block_ram_calibration_brick(n_samples; kwargs...)

Sample the standalone BRICK Turing model (`brick_model`) with a per-block Robust Adaptive
Metropolis sampler (see module docstring above): glaciers, thermal expansion, Antarctic ice sheet,
Greenland ice sheet, and the statistical auto-covariance parameters each get their own
adaptively-tuned multivariate proposal; `R` (cross-covariance) is a simple LKJ(4,η) prior-redraw
Metropolis step. Mirrors `run_per_block_ram_calibration` (the SNEASY-BRICK version).

Function Arguments:

    n_samples                = Number of outer sweeps (one RAM step per block, per sweep) to draw.
    output_dir                = Directory the chain (and a timing summary) are saved to.
    calibration_start_year   = First year to run the model (default 1850).
    calibration_end_year     = Final year of the calibration (default 2017).
    joint_antarctic_prior    = TRUE/FALSE: joint normal (TRUE) vs. marginal kernel density (FALSE, default)
                                Antarctic ice sheet prior.
    eta                       = Fixed `LKJ(4, η)` concentration for the `R` prior. `nothing` (default) chooses `η`
                                from the empirically estimated cross-stream correlation.
    reference_posterior_path  = Optional path to a CSV of posterior parameter draws (e.g. a previous
                                `run_calibration` subsample, such as `parameters_subsample_brick.csv`)
                                used to build each block's initial proposal covariance. Defaults to `nothing`,
                                which instead seeds every block from the package's own bundled
                                `initial_proposal_covariance_matrix_brick_04Jun2026.csv`. Either way this
                                is only a warm start -- RAM adapts away a poorly-scaled starting covariance
                                over the run.
    calibration_data_dir      = Directory containing the package's calibration data files. Defaults to the
                                bundled `data/calibration_data` directory.
    seed                      = RNG seed for reproducibility.

Returns `(chain_df, elapsed_seconds, acceptance_rates)`, where `chain_df` is a `DataFrame` with one
column per (scalar or vector-element) parameter and one row per sweep, and `acceptance_rates` is a
`Dict{Symbol,Float64}` of each block's overall acceptance rate.
"""
function run_per_block_ram_calibration_brick(n_samples::Int;
                                              output_dir::String,
                                              calibration_start_year::Int = 1850,
                                              calibration_end_year::Int = 2017,
                                              joint_antarctic_prior::Bool = false,
                                              eta::Union{Real, Nothing} = nothing,
                                              reference_posterior_path::Union{String, Nothing} = nothing,
                                              calibration_data_dir::Union{String, Nothing} = nothing,
                                              seed::Int = 2024)

    mkpath(output_dir)

    if isnothing(calibration_data_dir)
        calibration_data_dir = joinpath(@__DIR__, "..", "..", "data", "calibration_data")
    end

    run_mymodel! = construct_run_brick(calibration_start_year, calibration_end_year)
    model, empirical_R, eta = construct_brick_turing_model(run_mymodel!;
        model_start_year = calibration_start_year, calibration_end_year = calibration_end_year,
        joint_antarctic_prior = joint_antarctic_prior, eta = eta,
        calibration_data_dir = calibration_data_dir)
    println("Using η = $eta for the R ~ LKJ(4, η) prior.")

    #---------------------------------------------------------------------------
    # Seed the chain from the package's known-good calibration starting values (raw prior draws are
    # essentially always numerically invalid for this model).
    #---------------------------------------------------------------------------
    initvals = DataFrame(load(joinpath(calibration_data_dir, "calibration_initial_values_brick_04Jun2026.csv"))).parameter_values
    state = Dict{Symbol, Any}(zip(BRICK_SCALAR_PARAM_NAMES, initvals[1:20]))
    state[:antarctic_params] = Vector{Float64}(initvals[21:35])
    # Start `R` at the identity (no cross-stream correlation) rather than the empirical estimate --
    # see `run_per_block_ram_calibration`'s identical rationale: `empirical_R` combined with the
    # calibration_initial_values σ/ρ can make Σ_total non-PD at this exact point, giving -Inf at
    # sweep 1 and forcing every *other* block to reject until R happens to redraw into a compatible
    # state. Identity reliably gives a finite starting log-density.
    state[:R] = Matrix{Float64}(I, 4, 4)

    #---------------------------------------------------------------------------
    # Initial per-block proposal covariances (a warm start only -- RAM adapts away a poorly-scaled
    # guess over the run). Default: slice the package's own bundled covariance matrix, which shares
    # row/column order and units (including log-space `antarctic_precip0`) with
    # `calibration_initial_values_brick_04Jun2026.csv`. If the caller supplies
    # `reference_posterior_path` (e.g. a previous calibration's posterior subsample), use its
    # empirical covariance instead; that file stores `antarctic_precip0` in raw physical units, so
    # it's log-transformed first to match the model's internal log-space convention.
    #---------------------------------------------------------------------------
    if reference_posterior_path === nothing
        full_cov = Matrix(DataFrame(load(joinpath(calibration_data_dir, "initial_proposal_covariance_matrix_brick_04Jun2026.csv"))))
    else
        ram_df = DataFrame(load(reference_posterior_path))
        ram_df.antarctic_precip0 = log.(ram_df.antarctic_precip0)
        full_cov = Matrix(cov(Matrix(Float64.(ram_df[!, collect(BRICK_EXTERNAL_PARAM_NAMES)]))))
    end
    initial_cov(idx) = Symmetric((2.38^2 / length(idx)) .* full_cov[idx, idx])

    blocks = NamedTuple(
        block_name => BlockState(names, length(flatten_block(state, names)), cholesky(initial_cov(block_param_indices_brick(names))).L)
        for (block_name, names) in pairs(BRICK_GIBBS_BLOCKS)
    )

    #---------------------------------------------------------------------------
    # Log-density evaluator: rebuilds a VarInfo from the full named-parameter Dict every call (so
    # this pays the same per-evaluation cost -- one BRICK simulation -- as a Gibbs/MH approach; no
    # efficiency lost by not using Turing's own Gibbs machinery here).
    #---------------------------------------------------------------------------
    function logdensity(s::Dict{Symbol, Any})
        nt = NamedTuple(s)
        vi = DynamicPPL.VarInfo(model, DynamicPPL.InitFromParams(nt, nothing))
        lp = DynamicPPL.getlogp(vi)
        return lp.logprior + lp.logjac + lp.loglikelihood
    end

    #---------------------------------------------------------------------------
    # R: a simple (non-adaptive) LKJ(4,η) prior-redraw Metropolis step.
    #---------------------------------------------------------------------------
    r_prior = LKJ(4, eta)
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
    scalar_cols = Dict{Symbol, Vector{Float64}}(name => Vector{Float64}(undef, n_samples) for name in BRICK_SCALAR_PARAM_NAMES)
    antarctic_cols = [Vector{Float64}(undef, n_samples) for _ in 1:15]
    r_cols = [Vector{Float64}(undef, n_samples) for _ in 1:16]
    logjoint_col = Vector{Float64}(undef, n_samples)

    println("Sampling $n_samples sweeps with per-block RAM (glaciers / thermal / Antarctic / Greenland / auto-covariance / cross-covariance) for standalone BRICK...")
    elapsed = @elapsed for iter in 1:n_samples
        for (block_name, block) in pairs(blocks)
            accepted, current_lp = ram_block_step!(logdensity, current_lp, state, block; iter = iter)
            n_accept[block_name] += accepted
        end
        accepted_r, current_lp = r_step!(current_lp, state)
        n_accept[:R] += accepted_r

        for name in BRICK_SCALAR_PARAM_NAMES
            scalar_cols[name][iter] = state[name]
        end
        for i in 1:15
            antarctic_cols[i][iter] = state[:antarctic_params][i]
        end
        for i in 1:4, j in 1:4
            r_cols[(i-1)*4 + j][iter] = state[:R][i, j]
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
    for name in BRICK_SCALAR_PARAM_NAMES
        df[!, string(name)] = scalar_cols[name]
    end
    for i in 1:15
        df[!, "antarctic_params[$i]"] = antarctic_cols[i]
    end
    for i in 1:4, j in 1:4
        df[!, "R_$(i)_$(j)"] = r_cols[(i-1)*4 + j]
    end
    df[!, "logjoint"] = logjoint_col

    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    basename_ = "brick_ram_chain_n$(n_samples)_$(timestamp)"

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
