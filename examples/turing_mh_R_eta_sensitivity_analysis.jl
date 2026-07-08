using MimiBRICK
using Turing
using DataFrames
using CSVFiles
using Statistics
using Serialization

include(joinpath(@__DIR__, "turing_mh_calibration_benchmark.jl")) # brings in `run_mh_calibration` (its own auto-run guard means this doesn't trigger a 10,000-sample run)

##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
## Sensitivity analysis: how much does the choice of `η` -- the `LKJ(7, η)` concentration
## parameter for the (uncertain, sampled) cross-stream correlation matrix `R` -- shift the
## posterior for the physical SNEASY and BRICK parameters?
##
## `η` controls how strongly `R` is pulled toward the identity (independence) vs. allowed to take
## on more extreme correlations: larger `η` concentrates more mass near independence, `η=1` is
## uniform over all valid correlation matrices (no information about typical correlation
## strength), `η<1` favors more extreme correlations. `construct_sneasybrick_turing_model`'s
## default chooses `η` to be informed by the empirically observed cross-stream correlation
## strength (see `empirical_cross_correlation`/`choose_lkj_eta`) rather than defaulting to the
## uninformative `η=1`. This script compares that default against a spread of alternative `η`
## values to see whether the resulting BRICK/SNEASY parameter posteriors are actually sensitive to
## that choice, or whether it doesn't matter much in practice.
##
## Caveat: this model's MH sampler mixes poorly (see the benchmark script) -- ESS on the order of
## 1-12 out of 500 samples for most parameters. With that little effective information per chain,
## small differences between `η` specifications below may just be MCMC noise, not a real
## sensitivity to `η`. Read the comparison as a coarse, directional check, not a precise one; a
## much longer per-specification chain (and/or better step-size tuning) would be needed to say
## more confidently whether an observed difference is real.
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------

"""
    run_eta_sensitivity_analysis(n_samples::Int=200; etas=nothing, output_dir=...)

Run the SNEASY-BRICK MH calibration once per `η` value in `etas`, and compare the resulting
posterior mean/std for every SNEASY and BRICK (including the 15 Antarctic ice sheet) parameter
across specifications.

Function Arguments:

    n_samples  = Number of MH samples to draw per `η` specification (kept modest by default --
                 see the caveat above about mixing quality).
    etas       = Vector of `η` values to compare. Defaults to `[0.5, 1.0, η_data_informed, 3.0, 10.0]`,
                 where `η_data_informed` is `construct_sneasybrick_turing_model`'s own default
                 (chosen from the empirical cross-stream correlation).
    output_dir = Directory chains and the summary table are saved to.

Returns `(summary_df, chains_by_eta)`.
"""
function run_eta_sensitivity_analysis(n_samples::Int=200;
                                       etas::Union{Vector{<:Real}, Nothing} = nothing,
                                       output_dir::String = joinpath(@__DIR__, "results", "eta_sensitivity"))

    mkpath(output_dir)

    if isnothing(etas)
        println("Determining the data-informed default η...")
        _, _, eta_data_informed = MimiBRICK.construct_sneasybrick_turing_model(
            MimiBRICK.construct_run_sneasybrick(1850, 2017); model_start_year = 1850, calibration_end_year = 2017)
        etas = sort(unique([0.5, 1.0, round(eta_data_informed, digits = 3), 3.0, 10.0]))
        println("η specifications to compare: $etas  (data-informed default: $(round(eta_data_informed, digits=3)))")
    end

    sneasy_vars = ("CO₂_0", "N₂O_0", "temperature_0", "ocean_heat_0", "Q10", "CO₂_fertilization",
                   "CO₂_diffusivity", "heat_diffusivity", "rf_scale_aerosol", "ECS")
    brick_vars  = ("thermal_s₀", "greenland_v₀", "glaciers_v₀", "glaciers_s₀", "antarctic_s₀", "thermal_α",
                   "greenland_a", "greenland_b", "greenland_α", "greenland_β", "glaciers_β₀", "glaciers_n",
                   ("antarctic_params[$i]" for i in 1:15)...)

    summary_rows = NamedTuple[]
    chains_by_eta = Dict{Float64, Any}()

    for eta in etas
        println("="^70)
        println("Running calibration with η = $eta ($n_samples samples)...")
        chain, elapsed, _ = run_mh_calibration(n_samples; eta = eta, output_dir = output_dir, seed = 2024)
        chains_by_eta[eta] = chain
        println("η = $eta done in $(round(elapsed, digits=1))s.")

        for (block, names) in (("SNEASY", sneasy_vars), ("BRICK", brick_vars))
            for name in names
                samples = vec(Array(chain[:, name, :]))
                push!(summary_rows, (eta = eta, block = block, parameter = name,
                                      mean = mean(samples), std = std(samples), n_unique = length(unique(samples))))
            end
        end
    end

    summary_df = DataFrame(summary_rows)
    save(joinpath(output_dir, "eta_sensitivity_summary.csv"), summary_df)
    println()
    println("Saved per-(η, parameter) posterior mean/std/n_unique to eta_sensitivity_summary.csv")

    #---------------------------------------------------------------------------
    # Pivot for easy visual comparison: one row per parameter, one column per η.
    #---------------------------------------------------------------------------
    println()
    println("Posterior MEAN by η (SNEASY then BRICK parameters):")
    mean_pivot = unstack(summary_df, [:block, :parameter], :eta, :mean)
    show(mean_pivot, allrows = true, allcols = true)
    println()
    save(joinpath(output_dir, "eta_sensitivity_mean_pivot.csv"), mean_pivot)

    println()
    println()
    println("Posterior STD by η:")
    std_pivot = unstack(summary_df, [:block, :parameter], :eta, :std)
    show(std_pivot, allrows = true, allcols = true)
    println()
    save(joinpath(output_dir, "eta_sensitivity_std_pivot.csv"), std_pivot)

    return summary_df, chains_by_eta
end

#-------------------------------------------------------------------------------
# Default: 200 samples per η specification. Override by passing a sample count on the command
# line (`julia --project=. examples/turing_mh_R_eta_sensitivity_analysis.jl 500`) or by calling
# `run_eta_sensitivity_analysis(N; etas=[...])` directly from the REPL.
#-------------------------------------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    n_samples = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 200
    run_eta_sensitivity_analysis(n_samples)
end
