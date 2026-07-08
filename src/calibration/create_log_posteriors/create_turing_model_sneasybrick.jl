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
# SNEASY-BRICK model defined in `create_log_posterior_sneasybrick.jl`. Rather
# than building a bare log-posterior function for use with the package's
# RobustAdaptiveMetropolisSampler-based `run_calibration`, `construct_sneasybrick_turing_model`
# returns a Turing `Model` that can be sampled directly with any Turing sampler
# (e.g. `sample(model, NUTS(), 1000)`).
#-------------------------------------------------------------------------------

"""
    construct_sneasybrick_antarctic_priors(joint_antarctic_prior::Bool; calibration_data_dir::Union{String, Nothing} = nothing)

Load the paleo-calibrated Antarctic ice sheet parameters and build the pieces needed to score
the (informative) Antarctic ice sheet prior: parameter bounds and either a fitted joint
multivariate normal distribution or 15 marginal (truncated kernel density) distributions.

Description: This mirrors the Antarctic ice sheet prior construction from `construct_sneasybrick_log_prior`,
             factored out so it can be run once (outside of the Turing model body) and passed into
             `construct_sneasybrick_turing_model`.
"""
function construct_sneasybrick_antarctic_priors(joint_antarctic_prior::Bool; calibration_data_dir::Union{String, Nothing} = nothing)

    if isnothing(calibration_data_dir)
        calibration_data_dir = joinpath(@__DIR__, "..", "..", "..", "data", "calibration_data")
    end

    # Load required data to create Antarctic ice sheet informative priors (posterior parameters from previous calibration to paleo data).
    # Note: this excludes the Antarctic variance term because the model uses an AR(1) model for the recent instrumental observations.
    #       From original BRICK Fortran/R code: "var.dais was fit to paleo data-model mismatch, not representative of the current era."
    antarctic_paleo_file   = joinpath(calibration_data_dir, "DAISfastdyn_calibratedParameters_gamma_29Jan2017.nc")
    antarctic_paleo_params = convert(Array{Float64,2}, ncread(antarctic_paleo_file, "DAIS_parameters"))'[:,1:15]

    # Calculate upper and lower bounds for Antarctic ice sheet parameters (min/max values from paleo calibration).
    antarctic_lower_bound = vec(minimum(antarctic_paleo_params, dims=1))
    antarctic_upper_bound = vec(maximum(antarctic_paleo_params, dims=1))

    # Handle precip0 differently, since now sampling from log(P0)
    antarctic_lower_bound[7] = log(antarctic_lower_bound[7])
    antarctic_upper_bound[7] = log(antarctic_upper_bound[7])

    if joint_antarctic_prior
        # Fit a multivariate normal to the data (do not use variance estimate for paleo data, separately estimate AR(1) parameters for recent observations instead).
        antarctic_joint_prior    = fit(MvNormal, antarctic_paleo_params')
        antarctic_marginal_priors = nothing
    else
        # If not using joint prior, create marginal kernel density estimates for each Antarctic ice sheet model parameter.
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
    empirical_cross_correlation(f_run_model!, reference_params, n, calibration_data,
                                 indices_lawdome_co2_data, indices_maunaloa_co2_data, n_lawdome_co2, indices_co2_data,
                                 indices_temperature_data, indices_oceanheat_data, indices_glaciers_data,
                                 indices_greenland_data, indices_antarctic_data, indices_gmsl_data)

Estimate the fixed 7x7 cross-stream correlation matrix `R` used in `sneasybrick_model`'s joint
likelihood from data, rather than treating it as an uncertain/sampled parameter.

Description: Runs SNEASY-BRICK once at `reference_params` (a plausible, already-calibrated
             parameter vector -- see `construct_sneasybrick_turing_model`, which defaults to the
             package's `calibration_initial_values_sneasybrick_04Jun2026.csv`), computes each of
             the seven correlated streams' data-model residuals (temperature, ocean heat, CO₂,
             glaciers, greenland, antarctic, gmsl), and estimates their pairwise Pearson
             correlations at the years each pair of streams observes in common. Filling in
             pairwise correlations independently like this does not generally give a valid
             (positive semi-definite) correlation matrix once more than 2 streams are involved, so
             the result is projected onto the nearest valid correlation matrix (eigenvalue
             clipping) before being returned.
"""
function empirical_cross_correlation(f_run_model!, reference_params::Vector{Float64}, n::Int, calibration_data::DataFrame,
                                      indices_lawdome_co2_data::Vector{Int}, indices_maunaloa_co2_data::Vector{Int}, n_lawdome_co2::Int,
                                      indices_co2_data::Vector{Int}, indices_temperature_data::Vector{Int}, indices_oceanheat_data::Vector{Int},
                                      indices_glaciers_data::Vector{Int}, indices_greenland_data::Vector{Int},
                                      indices_antarctic_data::Vector{Int}, indices_gmsl_data::Vector{Int})

    modeled_CO₂               = Vector{Float64}(undef, n)
    modeled_oceanCO₂_flux     = Vector{Float64}(undef, n)
    modeled_temperature       = Vector{Float64}(undef, n)
    modeled_ocean_heat        = Vector{Float64}(undef, n)
    modeled_glaciers          = Vector{Float64}(undef, n)
    modeled_greenland         = Vector{Float64}(undef, n)
    modeled_antarctic         = Vector{Float64}(undef, n)
    modeled_thermal_expansion = Vector{Float64}(undef, n)
    modeled_gmsl              = Vector{Float64}(undef, n)
    f_run_model!(reference_params, modeled_CO₂, modeled_oceanCO₂_flux, modeled_temperature, modeled_ocean_heat,
                 modeled_glaciers, modeled_greenland, modeled_antarctic, modeled_thermal_expansion, modeled_gmsl)

    temperature_residual = Vector{Float64}(calibration_data[indices_temperature_data, :hadcrut_temperature_obs]) .- modeled_temperature[indices_temperature_data]
    ocean_heat_residual   = Vector{Float64}(calibration_data[indices_oceanheat_data, :ocean_heat_obs]) .- modeled_ocean_heat[indices_oceanheat_data]

    co2_residual = Vector{Float64}(undef, length(indices_co2_data))
    for (i, index) in enumerate(indices_lawdome_co2_data)
        co2_residual[i] = calibration_data[index, :lawdome_co2_obs] - mean(modeled_CO₂[index .+ (-4:3)])
    end
    for (i, index) in enumerate(indices_maunaloa_co2_data)
        co2_residual[i+n_lawdome_co2] = calibration_data[index, :maunaloa_co2_obs] - modeled_CO₂[index]
    end

    glaciers_residual   = Vector{Float64}(calibration_data[indices_glaciers_data, :glaciers_obs]) .- modeled_glaciers[indices_glaciers_data]
    greenland_residual  = Vector{Float64}(calibration_data[indices_greenland_data, :merged_greenland_obs]) .- modeled_greenland[indices_greenland_data]
    antarctic_residual  = Vector{Float64}(calibration_data[indices_antarctic_data, :antarctic_imbie_obs]) .- modeled_antarctic[indices_antarctic_data]
    gmsl_residual       = Vector{Float64}(calibration_data[indices_gmsl_data, :gmsl_obs]) .- modeled_gmsl[indices_gmsl_data]

    residual_streams    = (temperature_residual, ocean_heat_residual, co2_residual, glaciers_residual, greenland_residual, antarctic_residual, gmsl_residual)
    stream_year_indices = (indices_temperature_data, indices_oceanheat_data, indices_co2_data, indices_glaciers_data, indices_greenland_data, indices_antarctic_data, indices_gmsl_data)

    R = Matrix{Float64}(I, 7, 7)
    for a in 1:6, b in (a+1):7
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
    choose_lkj_eta(target_mean_abs_corr::Real; d::Int=7, n_mc::Int=5000, eta_grid=exp.(range(log(0.1), log(30), length=60)))

Choose the `LKJ(d, η)` concentration parameter `η` whose Monte Carlo-estimated average absolute
off-diagonal correlation is closest to `target_mean_abs_corr`.

Description: `LKJ(d, η)` has no built-in way to center its prior mass on an arbitrary target
             correlation matrix (larger `η` concentrates toward the identity/independence, `η=1`
             is uniform over all valid correlation matrices, `η<1` favors more extreme
             correlations) -- there's no direct "target matrix" parameter to set. This instead
             calibrates `η` so that random draws from `LKJ(d, η)` have, on average, about as much
             off-diagonal correlation as was empirically observed, which is a coarser (moment-
             matching), but tractable, way to make the *prior*'s typical correlation strength
             informed by the data rather than assuming no information about it (`η=1`).
"""
function choose_lkj_eta(target_mean_abs_corr::Real; d::Int=7, n_mc::Int=5000, eta_grid=exp.(range(log(0.1), log(30), length=60)))

    offdiag_idx = [(i, j) for i in 1:d for j in 1:d if i < j]

    best_eta, best_diff = eta_grid[1], Inf
    for η in eta_grid
        dist = LKJ(d, η)
        mean_abs_corr = mean(mean(abs(rand(dist)[i, j]) for (i, j) in offdiag_idx) for _ in 1:n_mc)
        diff = abs(mean_abs_corr - target_mean_abs_corr)
        if diff < best_diff
            best_eta, best_diff = η, diff
        end
    end

    return best_eta
end

"""
    sneasybrick_model(f_run_model!, n, calibration_data, indices_co2_data, indices_lawdome_co2_data,
                       indices_maunaloa_co2_data, n_lawdome_co2, indices_oceanco2_flux_data, indices_temperature_data,
                       indices_oceanheat_data, indices_glaciers_data, indices_greenland_data, indices_antarctic_data,
                       indices_gmsl_data, obs_thermal_trends, model_start_year, calibration_end_year,
                       antarctic_lower_bound, antarctic_upper_bound, joint_antarctic_prior, antarctic_joint_prior,
                       antarctic_marginal_priors, uniform_ECS)

Turing model for the SNEASY-BRICK calibration.

Description: Declares the same priors used by `construct_sneasybrick_log_prior` as Turing `~` statements
             (with the informative, non-Distributions.jl Antarctic ice sheet prior added manually via
             `Turing.@addlogprob!`), runs SNEASY-BRICK via `f_run_model!`, and scores all of the target
             outputs (temperature, ocean heat, CO₂, ocean CO₂ flux, glaciers, greenland, antarctic, thermal
             expansion trend, and GMSL) jointly under a single multivariate normal, rather than as separate
             independent data-set log-likelihood terms. Each data set keeps its own within-series AR(1)/CAR(1)/
             independent covariance structure, and an uncertain (sampled) cross-stream correlation matrix
             `R ~ LKJ(7, η)` additionally correlates the seven data sets that have their own fitted
             process-noise scale, with that cross-stream correlation decaying over the year gap between any
             pair of observations (not just contemporaneous ones) at a geometric mean of the two streams' own
             one-year-lag decay rates. `η` is a fixed hyperparameter (not sampled) -- see
             `construct_sneasybrick_turing_model` and `choose_lkj_eta`, which chooses a default `η` informed by
             `empirical_cross_correlation`'s data-derived estimate of the typical cross-stream correlation
             strength, rather than defaulting to an uninformative `η=1`.
"""
@model function sneasybrick_model(
    f_run_model!,
    n::Int,
    calibration_data::DataFrame,
    indices_co2_data::Vector{Int},
    indices_lawdome_co2_data::Vector{Int},
    indices_maunaloa_co2_data::Vector{Int},
    n_lawdome_co2::Int,
    indices_oceanco2_flux_data::Vector{Int},
    indices_temperature_data::Vector{Int},
    indices_oceanheat_data::Vector{Int},
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
    uniform_ECS::Bool,
    η::Real,
)

    #-----------------------------------------
    # Statistical Process Priors.
    #-----------------------------------------
    σ_temperature      ~ Uniform(1e-10, 0.2)  # SC-CH4 and Klaus papers.
    σ_ocean_heat       ~ Uniform(1e-10, 4)    # From SC-CH4 paper
    σ_glaciers         ~ Uniform(1e-10, 0.0015) # Based on BRICK code.
    σ_greenland        ~ Uniform(1e-10, 0.002)  # Based on BRICK code.
    σ_antarctic        ~ Uniform(1e-10, 0.063)  # Based on BRICK code.
    σ_gmsl             ~ Uniform(1e-10, 0.05)   # Just setting the same as prior_σ_gmsl_1900 value from old BRICK code.
    σ²_white_noise_CO₂ ~ Uniform(1e-10, 200)    # SC-CH4 paper.

    ρ_temperature ~ Uniform(-0.99, 0.99)
    ρ_ocean_heat  ~ Uniform(-0.99, 0.99)
    ρ_glaciers    ~ Uniform(-0.99, 0.99)
    ρ_greenland   ~ Uniform(-0.99, 0.99)
    ρ_antarctic   ~ Uniform(-0.99, 0.99)
    ρ_gmsl        ~ truncated(Normal(0.8, 0.25), -1.0, 1.0)
    α₀_CO₂        ~ Uniform(0.01, 11.5) # SC-CH4 paper.

    # Cross-data-set correlation matrix for the joint likelihood below. `η` (fixed, not sampled)
    # is chosen -- see `choose_lkj_eta` -- to be informed by the empirically estimated cross-stream
    # correlation strength (`empirical_cross_correlation`) rather than defaulting to η=1 (uniform
    # over all valid correlation matrices, i.e. no information about typical correlation strength).
    R ~ LKJ(7, η)

    #-----------------------------------------
    # Initial Condition Priors.
    #-----------------------------------------
    CO₂_0         ~ Uniform(275, 281)
    N₂O_0         ~ Uniform(264, 282)
    temperature_0 ~ truncated(Normal(), -0.3, 0.3)
    ocean_heat_0  ~ Uniform(-100, 0)
    thermal_s₀    ~ Uniform(-0.0484, 0.0484) # BRICK defaults. Initial sea level rise due to thermal expansion designated in 1850 (m SLE).
    greenland_v₀  ~ Uniform(7.16, 7.56)
    glaciers_v₀   ~ Uniform(0.31, 0.53)
    glaciers_s₀   ~ Uniform(-0.0536, 0.0791)
    antarctic_s₀  ~ Uniform(-0.04755, 0.05585) # Informed by prior BRICK runs.

    #-----------------------------------------
    # Carbon Cycle Priors.
    #-----------------------------------------
    Q10               ~ Uniform(1.0, 5)
    CO₂_fertilization ~ Uniform(0.0, 1)
    CO₂_diffusivity   ~ LogNormal(3.4, 0.43)

    #-----------------------------------------
    # Climate & Radiative Forcing Priors.
    #-----------------------------------------
    heat_diffusivity ~ LogNormal(1.1, 0.3)
    rf_scale_aerosol ~ TriangularDist(0.0, 3.0, 1.0)

    # Decide whether to use a uniform or paleo-informed ECS prior.
    ECS ~ uniform_ECS ? Uniform(0.0, 10.0) : truncated(Cauchy(3.0, 2.0), 0.0, 10.0)

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
        # `construct_sneasybrick_log_prior`, needed because Distributions.jl has no truncated MvNormal).
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
        σ_temperature, σ_ocean_heat, σ_glaciers, σ_greenland, σ_antarctic, σ_gmsl, σ²_white_noise_CO₂,
        ρ_temperature, ρ_ocean_heat, ρ_glaciers, ρ_greenland, ρ_antarctic, ρ_gmsl, α₀_CO₂,
        CO₂_0, N₂O_0, temperature_0, ocean_heat_0, thermal_s₀, greenland_v₀, glaciers_v₀, glaciers_s₀, antarctic_s₀,
        Q10, CO₂_fertilization, CO₂_diffusivity,
        heat_diffusivity, rf_scale_aerosol, ECS,
        thermal_α,
        greenland_a, greenland_b, greenland_α, greenland_β,
        glaciers_β₀, glaciers_n,
        antarctic_params...,
    ]
    T = eltype(p)

    #-----------------------------------------------------------------------
    # Run SNEASY+BRICK and calculate the data-model log-likelihood.
    # In case a parameter sample leads to non-physical model outcomes, reject rather than erroring out.
    #-----------------------------------------------------------------------
    modeled_CO₂               = Vector{T}(undef, n)
    modeled_oceanCO₂_flux     = Vector{T}(undef, n)
    modeled_temperature       = Vector{T}(undef, n)
    modeled_ocean_heat        = Vector{T}(undef, n)
    modeled_glaciers          = Vector{T}(undef, n)
    modeled_greenland         = Vector{T}(undef, n)
    modeled_antarctic         = Vector{T}(undef, n)
    modeled_thermal_expansion = Vector{T}(undef, n)
    modeled_gmsl              = Vector{T}(undef, n)

    try
        f_run_model!(p, modeled_CO₂, modeled_oceanCO₂_flux, modeled_temperature, modeled_ocean_heat,
                     modeled_glaciers, modeled_greenland, modeled_antarctic, modeled_thermal_expansion, modeled_gmsl)
    catch
        Turing.@addlogprob! -Inf
        return
    end

    #-----------------------------------------------------------------------------------------------------
    # Build a single joint multivariate normal log-likelihood spanning every target output. Each data set's
    # residuals keep their own (AR(1), CAR(1), or independent) covariance structure, but rather than summing
    # per-data-set log-likelihoods, all residuals/covariances are concatenated into one residual vector and
    # one block-diagonal covariance matrix so the whole observation set is scored under a single MvNormal.
    #
    # Everything below is wrapped in a single try/catch: besides the non-PD-covariance case (handled
    # explicitly), an out-of-range parameter combination reaching this point (e.g. |ρ| slightly outside
    # (-1,1) from a random-walk MH proposal that isn't support-aware) can throw a DomainError (sqrt of a
    # negative number in `process_std`/`decay_rate` below) — reject those draws instead of crashing the chain.
    #-----------------------------------------------------------------------------------------------------
    try

    #---------------------------------------------------------------------------
    # Global Surface Temperature (normalized to 1861-1880 mean).
    #---------------------------------------------------------------------------
    temperature_residual = Vector{Float64}(calibration_data[indices_temperature_data, :hadcrut_temperature_obs]) .- modeled_temperature[indices_temperature_data]
    Σ_temperature = ar1_covariance_matrix(length(temperature_residual), σ_temperature, ρ_temperature, calibration_data[indices_temperature_data, :hadcrut_temperature_sigma])

    #-----------------------------------------------------------------------
    # Ocean Heat Content.
    #-----------------------------------------------------------------------
    ocean_heat_residual = Vector{Float64}(calibration_data[indices_oceanheat_data, :ocean_heat_obs]) .- modeled_ocean_heat[indices_oceanheat_data]
    Σ_ocean_heat = ar1_covariance_matrix(length(ocean_heat_residual), σ_ocean_heat, ρ_ocean_heat, calibration_data[indices_oceanheat_data, :ocean_heat_sigma])

    #-----------------------------------------------------------------------
    # Atmospheric CO₂ Concentration.
    #-----------------------------------------------------------------------
    co2_residual = Vector{T}(undef, length(indices_co2_data))

    # Calculate CO₂ concentration (Law Dome) residuals (assuming 8 year model mean centered on year of ice core observation).
    for (i, index) in enumerate(indices_lawdome_co2_data)
        co2_residual[i] = calibration_data[index, :lawdome_co2_obs] - mean(modeled_CO₂[index .+ (-4:3)])
    end

    # Calculate CO₂ concentration (Mauna Loa) residuals.
    for (i, index) in enumerate(indices_maunaloa_co2_data)
        co2_residual[i+n_lawdome_co2] = calibration_data[index, :maunaloa_co2_obs] - modeled_CO₂[index]
    end

    Σ_co2 = car1_covariance_matrix(indices_co2_data, σ²_white_noise_CO₂, α₀_CO₂, calibration_data[indices_co2_data, :co2_combined_sigma])

    #-----------------------------------------------------------------------
    # Ocean Carbon Flux (independent observations).
    #-----------------------------------------------------------------------
    oceanco2_flux_residual = Vector{Float64}(calibration_data[indices_oceanco2_flux_data, :oceanco2_flux_obs]) .- modeled_oceanCO₂_flux[indices_oceanco2_flux_data]
    Σ_oceanco2_flux = Diagonal(Float64.(calibration_data[indices_oceanco2_flux_data, :oceanco2_flux_sigma]) .^ 2)

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
    # covariance matrix (each data set keeps its own AR(1)/CAR(1)/independent temporal
    # covariance structure on the block diagonal).
    #-----------------------------------------------------------------------
    total_residual = vcat(temperature_residual, ocean_heat_residual, co2_residual, oceanco2_flux_residual,
                           glaciers_residual, greenland_residual, antarctic_residual, thermal_trend_residual, gmsl_residual)

    Σ_total = Matrix(cat(Σ_temperature, Σ_ocean_heat, Σ_co2, Σ_oceanco2_flux,
                          Σ_glaciers, Σ_greenland, Σ_antarctic, Σ_thermal_trend, Σ_gmsl; dims=(1,2)))

    #-----------------------------------------------------------------------
    # Fill in cross-data-set covariance between *every* pair of positions across the seven
    # correlated streams (not just years they both happen to observe), using `R` for the
    # lag-0 correlation strength and each stream's own stationary process standard deviation.
    # The correlation between two streams' observations decays with the gap in years between
    # them — analogous to each stream's own within-series AR(1)/CAR(1) kernel — at a geometric-
    # mean of the two streams' own one-year-lag decay rates (|ρ| for the AR(1) streams, exp(-α₀)
    # for CO₂'s CAR(1) kernel). This is a simplifying, non-data-derived choice for the cross-stream
    # lag structure (no cross-stream lag data exists to fit it to), not a first-principles result.
    #-----------------------------------------------------------------------
    process_std = (
        σ_temperature / sqrt(1 - ρ_temperature^2),
        σ_ocean_heat  / sqrt(1 - ρ_ocean_heat^2),
        sqrt(σ²_white_noise_CO₂ / (2 * α₀_CO₂)),
        σ_glaciers    / sqrt(1 - ρ_glaciers^2),
        σ_greenland   / sqrt(1 - ρ_greenland^2),
        σ_antarctic   / sqrt(1 - ρ_antarctic^2),
        σ_gmsl        / sqrt(1 - ρ_gmsl^2),
    )

    decay_rate = (abs(ρ_temperature), abs(ρ_ocean_heat), exp(-α₀_CO₂), abs(ρ_glaciers), abs(ρ_greenland), abs(ρ_antarctic), abs(ρ_gmsl))

    stream_year_indices = (indices_temperature_data, indices_oceanheat_data, indices_co2_data,
                            indices_glaciers_data, indices_greenland_data, indices_antarctic_data, indices_gmsl_data)

    block_lengths = (length(temperature_residual), length(ocean_heat_residual), length(co2_residual), length(oceanco2_flux_residual),
                      length(glaciers_residual), length(greenland_residual), length(antarctic_residual), length(thermal_trend_residual), length(gmsl_residual))
    block_offsets = Tuple(cumsum((0, block_lengths[1:end-1]...)))
    stream_block   = (1, 2, 3, 5, 6, 7, 9) # each correlated stream's position among the 9 concatenated blocks

    for a in 1:6, b in (a+1):7
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
    # Score every target output under a single joint MvNormal. The ad hoc cross-block insertion
    # isn't guaranteed positive-definite for every (R, σ, ρ) draw -- replicating cross-stream
    # correlations out to every pair of residual positions does not itself guarantee the
    # resulting (much larger) joint covariance matrix stays positive-definite, even when the
    # small 7x7 `R` it's built from does. Rather than repair a failed draw (tried previously: an
    # eigenvalue-clipping repair works, but its eigendecomposition of the full ~400+-dimensional
    # `Σ_total` was the dominant per-sample cost whenever `R` had any non-trivial correlation,
    # ~15-65x slower than rejecting outright), just reject it like the DomainError case above --
    # `R` is sampled here (see `R ~ LKJ(7, η)`), so the chain can still find valid, better-
    # conditioned `R` draws over time instead of being stuck with one problematic fixed `R`.
    #-----------------------------------------------------------------------
    Turing.@addlogprob! logpdf(MvNormal(Symmetric(Σ_total)), total_residual)

    catch
        Turing.@addlogprob! -Inf
    end
end

"""
    construct_sneasybrick_turing_model(f_run_model!; model_start_year::Int=1850, calibration_end_year::Int=2017,
                                        joint_antarctic_prior::Bool=false, uniform_ECS::Bool=false,
                                        calibration_data_dir::Union{String, Nothing} = nothing,
                                        reference_params::Union{Vector{Float64}, Nothing} = nothing,
                                        eta::Union{Real, Nothing} = nothing)

Build a Turing `Model` for the SNEASY-BRICK calibration.

Description: This loads the calibration data/observations and the Antarctic ice sheet paleo prior information
             (once, outside of the Turing model body), estimates the empirical cross-stream correlation
             (`empirical_cross_correlation`) and uses it to choose the `LKJ(7, η)` prior's concentration `η`
             (`choose_lkj_eta`, unless `eta` is given explicitly), and returns a `(; model, empirical_R, eta)`
             named tuple. `model` is a conditioned Turing `Model` (with `R` sampled, not fixed) that can be
             sampled directly, e.g. `sample(model, MH(), 1000)` (NUTS is not usable on this model -- see
             examples/turing_nuts_calibration_test.jl); `empirical_R` and `eta` are returned too since they're
             each useful on their own (e.g. as a good starting value for `R`, or for a sensitivity analysis
             across different `eta` choices).

Function Arguments:

    f_run_model            = A function that runs the specific climate model version and returns the output being calibrated to observations.
    model_start_year       = First year to run the model (not necessarily first year of the calibration if model initializes earlier).
    calibration_end_year   = The final year to run the model calibration (defaults to 2017).
    joint_antarctic_prior  = TRUE/FALSE check for whether to use a joint normal prior distribution (TRUE) or fitted
                              marginal kernel density estimates (FALSE) for the Antarctic ice sheet parameters.
    uniform_ECS             = TRUE/FALSE check for whether or not to use a uniform prior distribution for the equilibrium
                              climate sensitivity (true = use uniform).
    calibration_data_dir    = Data directory for calibration data. Defaults to package calibration data directory, changing this is not recommended.
    reference_params        = The 51-element parameter vector (same order as `p` in `sneasybrick_model`) used to compute
                              data-model residuals for `empirical_cross_correlation`. Defaults to the package's own
                              `calibration_initial_values_sneasybrick_04Jun2026.csv` (a previously calibrated point).
    eta                     = Fixed `LKJ(7, η)` concentration parameter. Defaults to `nothing`, which chooses `η`
                              from the empirical cross-correlation via `choose_lkj_eta`; pass a specific value
                              (e.g. `1.0` for an uninformative/uniform prior) to override that choice, e.g. for a
                              sensitivity analysis across different `η` specifications.
"""
function construct_sneasybrick_turing_model(f_run_model!; model_start_year::Int=1850, calibration_end_year::Int=2017, joint_antarctic_prior::Bool=false, uniform_ECS::Bool=false, calibration_data_dir::Union{String, Nothing} = nothing, reference_params::Union{Vector{Float64}, Nothing} = nothing, eta::Union{Real, Nothing} = nothing)

    # Create a vector of calibration years and calculate total number of years to run model.
    calibration_years = collect(model_start_year:calibration_end_year)
    n = length(calibration_years)

    # set calibration data directory if one was not provided ie. it is set as nothing
    if isnothing(calibration_data_dir)
        calibration_data_dir = joinpath(@__DIR__, "..", "..", "..", "data", "calibration_data")
    end

    if isnothing(reference_params)
        reference_params = Vector{Float64}(DataFrame(load(joinpath(calibration_data_dir, "calibration_initial_values_sneasybrick_04Jun2026.csv"))).parameter_values)
    end

    # Get Antarctic ice sheet paleo-informed prior pieces (bounds + joint or marginal density).
    antarctic_lower_bound, antarctic_upper_bound, antarctic_joint_prior, antarctic_marginal_priors =
        construct_sneasybrick_antarctic_priors(joint_antarctic_prior, calibration_data_dir=calibration_data_dir)

    # Load calibration data/observations.
    calibration_data, obs_antarctic_trends, obs_thermal_trends = MimiBRICK.load_calibration_data(model_start_year, calibration_end_year, last_sea_level_norm_year=1990, calibration_data_dir=calibration_data_dir)

    # Calculate indices for each year that has an observation in calibration data sets.
    indices_maunaloa_co2_data  = findall(x -> !ismissing(x), calibration_data.maunaloa_co2_obs)
    indices_lawdome_co2_data   = findall(x -> !ismissing(x), calibration_data.lawdome_co2_obs)
    indices_oceanco2_flux_data = findall(x -> !ismissing(x), calibration_data.oceanco2_flux_obs)
    indices_temperature_data   = findall(x -> !ismissing(x), calibration_data.hadcrut_temperature_obs)
    indices_oceanheat_data     = findall(x -> !ismissing(x), calibration_data.ocean_heat_obs)
    indices_glaciers_data      = findall(x -> !ismissing(x), calibration_data.glaciers_obs)
    indices_greenland_data     = findall(x -> !ismissing(x), calibration_data.merged_greenland_obs) # Use merged Greenland data.
    indices_antarctic_data     = findall(x -> !ismissing(x), calibration_data.antarctic_imbie_obs)
    indices_gmsl_data          = findall(x -> !ismissing(x), calibration_data.gmsl_obs)

    # Combine CO₂ indices from Law Dome and Mauna Loa observations.
    indices_co2_data = sort(vcat(indices_lawdome_co2_data, indices_maunaloa_co2_data))

    # Combine CO₂ measurement errors from Law Dome and Mauna Loa observations (just for convenience).
    calibration_data.co2_combined_sigma = calibration_data.lawdome_co2_sigma
    calibration_data.co2_combined_sigma[indices_maunaloa_co2_data] = calibration_data.maunaloa_co2_sigma[indices_maunaloa_co2_data]

    # Calculate number of ice core observations for CO₂ (used for indexing).
    n_lawdome_co2 = length(indices_lawdome_co2_data)

    # Estimate the empirical cross-stream correlation matrix from `reference_params`'s data-model
    # residuals, and use it to choose `η` (unless a specific value was given).
    empirical_R = empirical_cross_correlation(f_run_model!, reference_params, n, calibration_data,
        indices_lawdome_co2_data, indices_maunaloa_co2_data, n_lawdome_co2, indices_co2_data,
        indices_temperature_data, indices_oceanheat_data, indices_glaciers_data, indices_greenland_data,
        indices_antarctic_data, indices_gmsl_data)

    if isnothing(eta)
        offdiag_idx = [(i, j) for i in 1:7 for j in 1:7 if i < j]
        target_mean_abs_corr = mean(abs(empirical_R[i, j]) for (i, j) in offdiag_idx)
        eta = choose_lkj_eta(target_mean_abs_corr)
    end

    # Return the Turing model (conditioned on the calibration data), along with the empirical
    # correlation matrix and the chosen `η` (both useful on their own -- see docstring above).
    model = sneasybrick_model(
        f_run_model!,
        n,
        calibration_data,
        indices_co2_data,
        indices_lawdome_co2_data,
        indices_maunaloa_co2_data,
        n_lawdome_co2,
        indices_oceanco2_flux_data,
        indices_temperature_data,
        indices_oceanheat_data,
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
        uniform_ECS,
        eta,
    )

    return (; model, empirical_R, eta)
end

##------------------------------------------------------------------------------
## End
##------------------------------------------------------------------------------
