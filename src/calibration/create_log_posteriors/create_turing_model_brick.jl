using Missings
using DataFrames
using Distributions
using NetCDF
using KernelDensity
using LinearAlgebra
using Statistics
using CSVFiles
using Turing

#-------------------------------------------------------------------------------
# This file contains a Turing.jl reimplementation of the log-posterior for the
# standalone BRICK model defined in `create_log_posterior_brick.jl`, mirroring
# `create_turing_model_sneasybrick.jl`'s joint (cross-correlated) multivariate normal
# likelihood but for BRICK's four correlated sea-level-rise data sets (glaciers, Greenland,
# Antarctica, GMSL) rather than SNEASY-BRICK's seven. `construct_brick_turing_model` returns a
# Turing `Model` that can be sampled directly with any Turing sampler (e.g. `sample(model,
# NUTS(), 1000)`, though see `examples/turing_nuts_calibration_test.jl` for why NUTS doesn't
# actually work with Mimi-based models).
#-------------------------------------------------------------------------------

"""
    construct_brick_antarctic_priors(joint_antarctic_prior::Bool; calibration_data_dir::Union{String, Nothing} = nothing)

Load the paleo-calibrated Antarctic ice sheet parameters and build the pieces needed to score
the (informative) Antarctic ice sheet prior: parameter bounds and either a fitted joint
multivariate normal distribution or 15 marginal (truncated kernel density) distributions.

Description: This mirrors `construct_sneasybrick_antarctic_priors` (identical Antarctic ice sheet
             prior, independent of which climate model BRICK is coupled to) so this file stays
             self-contained, matching the parallel `construct_brick_log_prior`/
             `construct_sneasybrick_log_prior` structure in the non-Turing log-posterior files.
"""
function construct_brick_antarctic_priors(joint_antarctic_prior::Bool; calibration_data_dir::Union{String, Nothing} = nothing)

    if isnothing(calibration_data_dir)
        calibration_data_dir = joinpath(@__DIR__, "..", "..", "..", "data", "calibration_data")
    end

    antarctic_paleo_file   = joinpath(calibration_data_dir, "DAISfastdyn_calibratedParameters_gamma_29Jan2017.nc")
    antarctic_paleo_params = convert(Array{Float64,2}, ncread(antarctic_paleo_file, "DAIS_parameters"))'[:,1:15]

    antarctic_lower_bound = vec(minimum(antarctic_paleo_params, dims=1))
    antarctic_upper_bound = vec(maximum(antarctic_paleo_params, dims=1))

    # Handle precip0 differently, since now sampling from log(P0)
    antarctic_lower_bound[7] = log(antarctic_lower_bound[7])
    antarctic_upper_bound[7] = log(antarctic_upper_bound[7])

    if joint_antarctic_prior
        antarctic_joint_prior    = fit(MvNormal, antarctic_paleo_params')
        antarctic_marginal_priors = nothing
    else
        antarctic_joint_prior = nothing
        antarctic_marginal_priors = [
            truncated_kernel(antarctic_paleo_params[:,1],       antarctic_lower_bound[1],  antarctic_upper_bound[1]),  # anto_α
            truncated_kernel(antarctic_paleo_params[:,2],       antarctic_lower_bound[2],  antarctic_upper_bound[2]),  # anto_β
            truncated_kernel(antarctic_paleo_params[:,3],       antarctic_lower_bound[3],  antarctic_upper_bound[3]),  # γ
            truncated_kernel(antarctic_paleo_params[:,4],       antarctic_lower_bound[4],  antarctic_upper_bound[4]),  # α
            truncated_kernel(antarctic_paleo_params[:,5],       antarctic_lower_bound[5],  antarctic_upper_bound[5]),  # μ
            truncated_kernel(antarctic_paleo_params[:,6],       antarctic_lower_bound[6],  antarctic_upper_bound[6]),  # ν
            truncated_kernel(log.(antarctic_paleo_params[:,7]), antarctic_lower_bound[7],  antarctic_upper_bound[7]),  # precip₀ (log space)
            truncated_kernel(antarctic_paleo_params[:,8],       antarctic_lower_bound[8],  antarctic_upper_bound[8]),  # κ
            truncated_kernel(antarctic_paleo_params[:,9],       antarctic_lower_bound[9],  antarctic_upper_bound[9]),  # flow₀
            truncated_kernel(antarctic_paleo_params[:,10],      antarctic_lower_bound[10], antarctic_upper_bound[10]), # runoff_height₀
            truncated_kernel(antarctic_paleo_params[:,11],      antarctic_lower_bound[11], antarctic_upper_bound[11]), # c
            truncated_kernel(antarctic_paleo_params[:,12],      antarctic_lower_bound[12], antarctic_upper_bound[12]), # bedheight₀
            truncated_kernel(antarctic_paleo_params[:,13],      antarctic_lower_bound[13], antarctic_upper_bound[13]), # slope
            truncated_kernel(antarctic_paleo_params[:,14],      antarctic_lower_bound[14], antarctic_upper_bound[14]), # λ
            truncated_kernel(antarctic_paleo_params[:,15],      antarctic_lower_bound[15], antarctic_upper_bound[15]), # temp_threshold
        ]
    end

    return antarctic_lower_bound, antarctic_upper_bound, antarctic_joint_prior, antarctic_marginal_priors
end

"""
    brick_empirical_cross_correlation(f_run_model!, reference_params, n, calibration_data,
                                       indices_glaciers_data, indices_greenland_data,
                                       indices_antarctic_data, indices_gmsl_data)

Estimate the fixed 4x4 cross-stream correlation matrix `R` used in `brick_model`'s joint
likelihood from data, rather than treating it as an uncertain/sampled parameter.

Description: Mirrors `empirical_cross_correlation` from `create_turing_model_sneasybrick.jl`, but
             for BRICK's four correlated sea-level-rise streams (glaciers, Greenland, Antarctica,
             GMSL) instead of SNEASY-BRICK's seven -- BRICK itself doesn't simulate temperature,
             ocean heat, or CO₂ (those are exogenous forcing here, not calibration targets).
"""
function brick_empirical_cross_correlation(f_run_model!, reference_params::Vector{Float64}, n::Int, calibration_data::DataFrame,
                                            indices_glaciers_data::Vector{Int}, indices_greenland_data::Vector{Int},
                                            indices_antarctic_data::Vector{Int}, indices_gmsl_data::Vector{Int})

    modeled_glaciers          = Vector{Float64}(undef, n)
    modeled_greenland         = Vector{Float64}(undef, n)
    modeled_antarctic         = Vector{Float64}(undef, n)
    modeled_thermal_expansion = Vector{Float64}(undef, n)
    modeled_gmsl              = Vector{Float64}(undef, n)
    f_run_model!(reference_params, modeled_glaciers, modeled_greenland, modeled_antarctic, modeled_thermal_expansion, modeled_gmsl)

    glaciers_residual   = Vector{Float64}(calibration_data[indices_glaciers_data, :glaciers_obs]) .- modeled_glaciers[indices_glaciers_data]
    greenland_residual  = Vector{Float64}(calibration_data[indices_greenland_data, :merged_greenland_obs]) .- modeled_greenland[indices_greenland_data]
    antarctic_residual  = Vector{Float64}(calibration_data[indices_antarctic_data, :antarctic_imbie_obs]) .- modeled_antarctic[indices_antarctic_data]
    gmsl_residual       = Vector{Float64}(calibration_data[indices_gmsl_data, :gmsl_obs]) .- modeled_gmsl[indices_gmsl_data]

    residual_streams    = (glaciers_residual, greenland_residual, antarctic_residual, gmsl_residual)
    stream_year_indices = (indices_glaciers_data, indices_greenland_data, indices_antarctic_data, indices_gmsl_data)

    R = Matrix{Float64}(I, 4, 4)
    for a in 1:3, b in (a+1):4
        common_years = intersect(stream_year_indices[a], stream_year_indices[b])
        length(common_years) < 2 && continue # not enough overlap to estimate a correlation; leave at 0

        pos_a = Dict(idx => i for (i, idx) in enumerate(stream_year_indices[a]))
        pos_b = Dict(idx => i for (i, idx) in enumerate(stream_year_indices[b]))
        xa = [residual_streams[a][pos_a[yr]] for yr in common_years]
        xb = [residual_streams[b][pos_b[yr]] for yr in common_years]

        corr = cor(xa, xb)
        R[a, b] = corr
        R[b, a] = corr
    end

    # Project onto the nearest valid correlation matrix (eigenvalue clipping), since independently
    # estimated pairwise correlations aren't guaranteed to form a positive semi-definite matrix.
    R_sym               = Symmetric((R .+ R') ./ 2)
    eigenvalues, eigvecs = eigen(R_sym)
    R_psd                = eigvecs * Diagonal(max.(eigenvalues, 1e-6)) * eigvecs'
    d                    = sqrt.(diag(R_psd))
    return Matrix(Symmetric(R_psd ./ (d * d')))
end

"""
    brick_model(f_run_model!, n, calibration_data, indices_glaciers_data, indices_greenland_data,
                indices_antarctic_data, indices_gmsl_data, obs_thermal_trends, model_start_year,
                calibration_end_year, antarctic_lower_bound, antarctic_upper_bound,
                joint_antarctic_prior, antarctic_joint_prior, antarctic_marginal_priors, η)

Turing model for the standalone BRICK calibration.

Description: Declares the same priors used by `construct_brick_log_prior` as Turing `~` statements
             (with the informative, non-Distributions.jl Antarctic ice sheet prior added manually via
             `Turing.@addlogprob!`), runs BRICK via `f_run_model!`, and scores glaciers, Greenland,
             Antarctica, and GMSL jointly under a single multivariate normal, rather than as separate
             independent data-set log-likelihood terms (mirroring `sneasybrick_model`). The thermal
             expansion trend observations stay an independent Normal term, same as in
             `construct_brick_log_posterior` and in `sneasybrick_model` -- there's no year-indexed
             per-observation residual to correlate them against the other four streams (they're
             1971-2009/1993-2009 trend estimates, not annual observations). Each of the four
             correlated data sets keeps its own within-series AR(1) covariance structure, and an
             uncertain (sampled) cross-stream correlation matrix `R ~ LKJ(4, η)` additionally
             correlates them, with that cross-stream correlation decaying over the year gap between
             any pair of observations (not just contemporaneous ones) at a geometric mean of the two
             streams' own one-year-lag decay rates. `η` is a fixed hyperparameter (not sampled) --
             see `construct_brick_turing_model` and `choose_lkj_eta` (from
             `create_turing_model_sneasybrick.jl`, reused as-is since it isn't specific to any one
             model configuration), which chooses a default `η` informed by
             `brick_empirical_cross_correlation`'s data-derived estimate of the typical cross-stream
             correlation strength, rather than defaulting to an uninformative `η=1`.
"""
@model function brick_model(
    f_run_model!,
    n::Int,
    calibration_data::DataFrame,
    indices_glaciers_data::Vector{Int},
    indices_greenland_data::Vector{Int},
    indices_antarctic_data::Vector{Int},
    indices_gmsl_data::Vector{Int},
    obs_thermal_trends,
    model_start_year::Int,
    calibration_end_year::Int,
    antarctic_lower_bound::Vector{Float64},
    antarctic_upper_bound::Vector{Float64},
    joint_antarctic_prior::Bool,
    antarctic_joint_prior,
    antarctic_marginal_priors,
    η::Real,
)

    #-----------------------------------------
    # Statistical Process Priors.
    #-----------------------------------------
    σ_glaciers  ~ Uniform(1e-10, 0.0015) # Based on BRICK code.
    σ_greenland ~ Uniform(1e-10, 0.002)  # Based on BRICK code.
    σ_antarctic ~ Uniform(1e-10, 0.063)  # Based on BRICK code.
    σ_gmsl      ~ Uniform(1e-10, 0.05)   # Just setting the same as prior_σ_gmsl_1900 value from old BRICK code.

    ρ_glaciers  ~ Uniform(-0.99, 0.99)
    ρ_greenland ~ Uniform(-0.99, 0.99)
    ρ_antarctic ~ Uniform(-0.99, 0.99)
    ρ_gmsl      ~ truncated(Normal(0.8, 0.25), -1.0, 1.0)

    # Cross-data-set correlation matrix for the joint likelihood below. `η` (fixed, not sampled)
    # is chosen -- see `choose_lkj_eta` -- to be informed by the empirically estimated cross-stream
    # correlation strength (`brick_empirical_cross_correlation`) rather than defaulting to η=1
    # (uniform over all valid correlation matrices, i.e. no information about typical correlation strength).
    R ~ LKJ(4, η)

    #-----------------------------------------
    # Initial Condition Priors.
    #-----------------------------------------
    thermal_s₀   ~ Uniform(-0.0484, 0.0484) # BRICK defaults. Initial sea level rise due to thermal expansion designated in 1850 (m SLE).
    greenland_v₀ ~ Uniform(7.16, 7.56)
    glaciers_v₀  ~ Uniform(0.31, 0.53)
    glaciers_s₀  ~ Uniform(-0.0536, 0.0791)
    antarctic_s₀ ~ Uniform(-0.04755, 0.05585) # Informed by prior BRICK runs.

    #---------------------------------------------
    # Sea Level Rise From Thermal Expansion Priors
    #---------------------------------------------
    thermal_α ~ Uniform(0.05, 0.3) # Global ocean-averaged thermal expansion coefficient (kg m⁻³ °C⁻¹).

    #------------------------------------------------
    # Sea Level Rise from Greenland Ice Sheet Priors
    #------------------------------------------------
    greenland_a ~ Uniform(-4.0, -0.001)
    greenland_b ~ Uniform(5.888, 8.832)
    greenland_α ~ Uniform(0.0, 0.001)
    greenland_β ~ Uniform(0.0, 0.001)

    #------------------------------------------------------
    # Sea Level Rise from Glaciers & Small Ice Caps Priors
    #------------------------------------------------------
    glaciers_β₀ ~ Uniform(0.0, 0.041)
    glaciers_n  ~ Uniform(0.55, 1.0)

    #------------------------------------------------------------------------------------
    # Antarctic Ice Sheet Priors (informed by a previous calibration to paleo data).
    # Distributions.jl has no truncated multivariate normal / kernel density type, so these
    # are sampled under a flat/unconstrained proposal and scored manually via `@addlogprob!`.
    #------------------------------------------------------------------------------------
    if joint_antarctic_prior
        antarctic_params ~ antarctic_joint_prior

        # If sampled parameters are outside paleo range, reject (mirrors the -Inf check in
        # `construct_brick_log_prior`, needed because Distributions.jl has no truncated MvNormal).
        if any(antarctic_params .< antarctic_lower_bound) || any(antarctic_params .> antarctic_upper_bound)
            Turing.@addlogprob! -Inf
            return
        end
    else
        antarctic_params ~ arraydist(fill(Turing.Flat(), 15))
        # `antarctic_marginal_priors[i]` are `KernelDensity.InterpKDE` objects, which support `pdf` but not `logpdf`.
        Turing.@addlogprob! sum(log(pdf(antarctic_marginal_priors[i], antarctic_params[i])) for i in 1:15)
    end

    #------------------------------------------------------------------------------------
    # Assemble the full parameter vector in the order expected by `f_run_model!`. Left
    # untyped (rather than `Float64[...]`) so that under AD samplers (e.g. NUTS) this stays
    # a `Vector{<:ForwardDiff.Dual}` instead of erroring by truncating derivative information.
    #------------------------------------------------------------------------------------
    p = [
        σ_glaciers, σ_greenland, σ_antarctic, σ_gmsl,
        ρ_glaciers, ρ_greenland, ρ_antarctic, ρ_gmsl,
        thermal_s₀, greenland_v₀, glaciers_v₀, glaciers_s₀, antarctic_s₀,
        thermal_α,
        greenland_a, greenland_b, greenland_α, greenland_β,
        glaciers_β₀, glaciers_n,
        antarctic_params...,
    ]
    T = eltype(p)

    #-----------------------------------------------------------------------
    # Run BRICK and calculate the data-model log-likelihood.
    # In case a parameter sample leads to non-physical model outcomes, reject rather than erroring out.
    #-----------------------------------------------------------------------
    modeled_glaciers          = Vector{T}(undef, n)
    modeled_greenland         = Vector{T}(undef, n)
    modeled_antarctic         = Vector{T}(undef, n)
    modeled_thermal_expansion = Vector{T}(undef, n)
    modeled_gmsl              = Vector{T}(undef, n)

    try
        f_run_model!(p, modeled_glaciers, modeled_greenland, modeled_antarctic, modeled_thermal_expansion, modeled_gmsl)
    catch
        Turing.@addlogprob! -Inf
        return
    end

    #-----------------------------------------------------------------------------------------------------
    # Build a single joint multivariate normal log-likelihood spanning glaciers, Greenland, Antarctica,
    # the (independent) thermal expansion trend, and GMSL. Each of the four correlated data sets keeps its
    # own AR(1) covariance structure, but rather than summing per-data-set log-likelihoods, all residuals/
    # covariances are concatenated into one residual vector and one block-diagonal covariance matrix so the
    # whole observation set is scored under a single MvNormal.
    #
    # Everything below is wrapped in a single try/catch: besides the non-PD-covariance case (handled
    # implicitly by `MvNormal` throwing), an out-of-range parameter combination reaching this point (e.g.
    # |ρ| slightly outside (-1,1) from a random-walk MH proposal that isn't support-aware) can throw a
    # DomainError (sqrt of a negative number in `process_std` below) -- reject those draws instead of
    # crashing the chain.
    #-----------------------------------------------------------------------------------------------------
    try

    #-----------------------------------------------------------------------
    # Glaciers and Small Ice Caps.
    #-----------------------------------------------------------------------
    glaciers_residual = Vector{Float64}(calibration_data[indices_glaciers_data, :glaciers_obs]) .- modeled_glaciers[indices_glaciers_data]
    Σ_glaciers = ar1_covariance_matrix(length(glaciers_residual), σ_glaciers, ρ_glaciers, calibration_data[indices_glaciers_data, :glaciers_sigma])

    #-------------------------------------------------------------------------------
    # Greenland Ice Sheet Merged Data (normalized to 1992_2001 mean).
    #-------------------------------------------------------------------------------
    greenland_residual = Vector{Float64}(calibration_data[indices_greenland_data, :merged_greenland_obs]) .- modeled_greenland[indices_greenland_data]
    Σ_greenland = ar1_covariance_matrix(length(greenland_residual), σ_greenland, ρ_greenland, calibration_data[indices_greenland_data, :merged_greenland_sigma])

    #------------------------------------------------------------------------------------
    # AIS (Antarctic Ice Sheet) IMBIE Data (normalized to 1992_2001 mean).
    #------------------------------------------------------------------------------------
    antarctic_residual = Vector{Float64}(calibration_data[indices_antarctic_data, :antarctic_imbie_obs]) .- modeled_antarctic[indices_antarctic_data]
    Σ_antarctic = ar1_covariance_matrix(length(antarctic_residual), σ_antarctic, ρ_antarctic, calibration_data[indices_antarctic_data, :antarctic_imbie_sigma])

    #-----------------------------------------------------------------------
    # Thermal Expansion Trends for 1971-2009 & 1993-2009 (independent observations).
    # Add a check for whether or not the calibration period covers the thermal trend period
    # (some sensitivity tests will not), otherwise this data set contributes nothing.
    #-----------------------------------------------------------------------
    if calibration_end_year >= 2009
        modeled_thermal_trend = calculate_trends(modeled_thermal_expansion, obs_thermal_trends, model_start_year, calibration_end_year)
        thermal_trend_residual = modeled_thermal_trend .- Vector{Float64}(obs_thermal_trends.Trend)
        # σ for observed trends based on IPCC 90% trend window values.
        obs_trend_err = 0.5 .* (obs_thermal_trends.Upper_90_Percent .- obs_thermal_trends.Lower_90_Percent)
        Σ_thermal_trend = Diagonal(Float64.(obs_trend_err) .^ 2)
    else
        thermal_trend_residual = T[]
        Σ_thermal_trend = Diagonal(Float64[])
    end

    #---------------------------------------------------------------------------
    # Global Mean Sea Level Rise (normalized to 1961-1990 mean).
    #---------------------------------------------------------------------------
    gmsl_residual = Vector{Float64}(calibration_data[indices_gmsl_data, :gmsl_obs]) .- modeled_gmsl[indices_gmsl_data]
    Σ_gmsl = ar1_covariance_matrix(length(gmsl_residual), σ_gmsl, ρ_gmsl, calibration_data[indices_gmsl_data, :gmsl_sigma])

    #-----------------------------------------------------------------------
    # Concatenate every target output into one residual vector and one block-diagonal
    # covariance matrix (each correlated data set keeps its own AR(1) temporal covariance
    # structure on the block diagonal; the thermal expansion trend block stays independent).
    #-----------------------------------------------------------------------
    total_residual = vcat(glaciers_residual, greenland_residual, antarctic_residual, thermal_trend_residual, gmsl_residual)

    Σ_total = Matrix(cat(Σ_glaciers, Σ_greenland, Σ_antarctic, Σ_thermal_trend, Σ_gmsl; dims=(1,2)))

    #-----------------------------------------------------------------------
    # Fill in cross-data-set covariance between *every* pair of positions across the four
    # correlated streams (not just years they both happen to observe), using `R` for the
    # lag-0 correlation strength and each stream's own stationary process standard deviation.
    # The correlation between two streams' observations decays with the gap in years between
    # them -- analogous to each stream's own within-series AR(1) kernel -- at a geometric mean
    # of the two streams' own one-year-lag decay rates (|ρ|). This is a simplifying, non-data-
    # derived choice for the cross-stream lag structure (no cross-stream lag data exists to fit
    # it to), not a first-principles result -- mirrors `sneasybrick_model`'s treatment.
    #-----------------------------------------------------------------------
    process_std = (
        σ_glaciers  / sqrt(1 - ρ_glaciers^2),
        σ_greenland / sqrt(1 - ρ_greenland^2),
        σ_antarctic / sqrt(1 - ρ_antarctic^2),
        σ_gmsl      / sqrt(1 - ρ_gmsl^2),
    )

    decay_rate = (abs(ρ_glaciers), abs(ρ_greenland), abs(ρ_antarctic), abs(ρ_gmsl))

    stream_year_indices = (indices_glaciers_data, indices_greenland_data, indices_antarctic_data, indices_gmsl_data)

    block_lengths = (length(glaciers_residual), length(greenland_residual), length(antarctic_residual), length(thermal_trend_residual), length(gmsl_residual))
    block_offsets = Tuple(cumsum((0, block_lengths[1:end-1]...)))
    stream_block   = (1, 2, 3, 5) # each correlated stream's position among the 5 concatenated blocks (4 = thermal trend, uncorrelated)

    for a in 1:3, b in (a+1):4
        pos_a, pos_b = stream_year_indices[a], stream_year_indices[b]
        lag  = abs.(pos_a .- pos_b') # length(pos_a) x length(pos_b) matrix of year gaps
        ρ_ab = sqrt(decay_rate[a] * decay_rate[b])
        cross_block = (R[a, b] * process_std[a] * process_std[b]) .* (ρ_ab .^ lag)

        off_a, off_b = block_offsets[stream_block[a]], block_offsets[stream_block[b]]
        rows_a, rows_b = (off_a+1):(off_a+length(pos_a)), (off_b+1):(off_b+length(pos_b))
        Σ_total[rows_a, rows_b] = cross_block
        Σ_total[rows_b, rows_a] = cross_block'
    end

    #-----------------------------------------------------------------------
    # Score every target output under a single joint MvNormal. `R` is sampled here (see
    # `R ~ LKJ(4, η)`), so the chain can still find valid, better-conditioned `R` draws over
    # time instead of being stuck with one problematic fixed `R` -- see `sneasybrick_model`'s
    # docstring for why a failed draw is rejected outright rather than repaired.
    #-----------------------------------------------------------------------
    Turing.@addlogprob! logpdf(MvNormal(Symmetric(Σ_total)), total_residual)

    catch
        Turing.@addlogprob! -Inf
    end
end

"""
    construct_brick_turing_model(f_run_model!; model_start_year::Int=1850, calibration_end_year::Int=2017,
                                  joint_antarctic_prior::Bool=false,
                                  calibration_data_dir::Union{String, Nothing} = nothing,
                                  reference_params::Union{Vector{Float64}, Nothing} = nothing,
                                  eta::Union{Real, Nothing} = nothing)

Build a Turing `Model` for the standalone BRICK calibration.

Description: This loads the calibration data/observations and the Antarctic ice sheet paleo prior information
             (once, outside of the Turing model body), estimates the empirical cross-stream correlation
             (`brick_empirical_cross_correlation`) and uses it to choose the `LKJ(4, η)` prior's concentration `η`
             (`choose_lkj_eta`, unless `eta` is given explicitly), and returns a `(; model, empirical_R, eta)`
             named tuple. `model` is a conditioned Turing `Model` (with `R` sampled, not fixed) that can be
             sampled directly, e.g. `sample(model, MH(), 1000)`; `empirical_R` and `eta` are returned too since
             they're each useful on their own (e.g. as a good starting value for `R`, or for a sensitivity
             analysis across different `eta` choices) -- mirrors `construct_sneasybrick_turing_model`.

Function Arguments:

    f_run_model            = A function that runs BRICK and returns the output being calibrated to observations
                              (i.e. `construct_run_brick`'s `run_brick!`).
    model_start_year       = First year to run the model (not necessarily first year of the calibration if model initializes earlier).
    calibration_end_year   = The final year to run the model calibration (defaults to 2017).
    joint_antarctic_prior  = TRUE/FALSE check for whether to use a joint normal prior distribution (TRUE) or fitted
                              marginal kernel density estimates (FALSE) for the Antarctic ice sheet parameters.
    calibration_data_dir    = Data directory for calibration data. Defaults to package calibration data directory, changing this is not recommended.
    reference_params        = The 35-element parameter vector (same order as `p` in `brick_model`) used to compute
                              data-model residuals for `brick_empirical_cross_correlation`. Defaults to the package's own
                              `calibration_initial_values_brick_04Jun2026.csv` (a previously calibrated point).
    eta                     = Fixed `LKJ(4, η)` concentration parameter. Defaults to `nothing`, which chooses `η`
                              from the empirical cross-correlation via `choose_lkj_eta`; pass a specific value
                              (e.g. `1.0` for an uninformative/uniform prior) to override that choice, e.g. for a
                              sensitivity analysis across different `η` specifications.
"""
function construct_brick_turing_model(f_run_model!; model_start_year::Int=1850, calibration_end_year::Int=2017, joint_antarctic_prior::Bool=false, calibration_data_dir::Union{String, Nothing} = nothing, reference_params::Union{Vector{Float64}, Nothing} = nothing, eta::Union{Real, Nothing} = nothing)

    # Create a vector of calibration years and calculate total number of years to run model.
    calibration_years = collect(model_start_year:calibration_end_year)
    n = length(calibration_years)

    # set calibration data directory if one was not provided ie. it is set as nothing
    if isnothing(calibration_data_dir)
        calibration_data_dir = joinpath(@__DIR__, "..", "..", "..", "data", "calibration_data")
    end

    if isnothing(reference_params)
        reference_params = Vector{Float64}(DataFrame(load(joinpath(calibration_data_dir, "calibration_initial_values_brick_04Jun2026.csv"))).parameter_values)
    end

    # Get Antarctic ice sheet paleo-informed prior pieces (bounds + joint or marginal density).
    antarctic_lower_bound, antarctic_upper_bound, antarctic_joint_prior, antarctic_marginal_priors =
        construct_brick_antarctic_priors(joint_antarctic_prior, calibration_data_dir=calibration_data_dir)

    # Load calibration data/observations.
    calibration_data, obs_antarctic_trends, obs_thermal_trends = MimiBRICK.load_calibration_data(model_start_year, calibration_end_year, last_sea_level_norm_year=1990, calibration_data_dir=calibration_data_dir)

    # Calculate indices for each year that has an observation in calibration data sets.
    indices_glaciers_data  = findall(x -> !ismissing(x), calibration_data.glaciers_obs)
    indices_greenland_data = findall(x -> !ismissing(x), calibration_data.merged_greenland_obs) # Use merged Greenland data.
    indices_antarctic_data = findall(x -> !ismissing(x), calibration_data.antarctic_imbie_obs)
    indices_gmsl_data      = findall(x -> !ismissing(x), calibration_data.gmsl_obs)

    # Estimate the empirical cross-stream correlation matrix from `reference_params`'s data-model
    # residuals, and use it to choose `η` (unless a specific value was given).
    empirical_R = brick_empirical_cross_correlation(f_run_model!, reference_params, n, calibration_data,
        indices_glaciers_data, indices_greenland_data, indices_antarctic_data, indices_gmsl_data)

    if isnothing(eta)
        offdiag_idx = [(i, j) for i in 1:4 for j in 1:4 if i < j]
        target_mean_abs_corr = mean(abs(empirical_R[i, j]) for (i, j) in offdiag_idx)
        eta = choose_lkj_eta(target_mean_abs_corr; d=4)
    end

    # Return the Turing model (conditioned on the calibration data), along with the empirical
    # correlation matrix and the chosen `η` (both useful on their own -- see docstring above).
    model = brick_model(
        f_run_model!,
        n,
        calibration_data,
        indices_glaciers_data,
        indices_greenland_data,
        indices_antarctic_data,
        indices_gmsl_data,
        obs_thermal_trends,
        model_start_year,
        calibration_end_year,
        antarctic_lower_bound,
        antarctic_upper_bound,
        joint_antarctic_prior,
        antarctic_joint_prior,
        antarctic_marginal_priors,
        eta,
    )

    return (; model, empirical_R, eta)
end

##------------------------------------------------------------------------------
## End
##------------------------------------------------------------------------------
