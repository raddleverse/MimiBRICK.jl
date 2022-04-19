# Load required data to create Antarctic ice sheet informative priors (posterior parameters from previous calibration to paleo data).
# Note: this excludes the Antarctic variance term because the model uses an AR(1) model for the recent instrumental observations.
#       From original BRICK Fortran/R code: "var.dais was fit to paleo data-model mismatch, not representative of the current era."
antarctic_paleo_file   = joinpath(@__DIR__, "..","data", "calibration_data", "DAISfastdyn_calibratedParameters_gamma_29Jan2017.nc")
antarctic_paleo_params = convert(Array{Float64,2}, ncread(antarctic_paleo_file, "DAIS_parameters"))'[:,1:15]

# #-------------------------------------------------------------------------------------------------------
# #-------------------------------------------------------------------------------------------------------
# # This file contains functions used to calculate the log-posterior for the SNEASY-BRICK model.
# #-------------------------------------------------------------------------------------------------------
# #-------------------------------------------------------------------------------------------------------


#######################################################################################################################
# CALCULATE TOTAL (LOG) PRIOR PROBABILITY.
#######################################################################################################################
# Description: This creates a function that will calculate the total (log) prior probability of the uncertain model,
#              initial condition, and statistical process parameters specific to the SNEASY-BRICK model. It uses
#              non-uniform priors for the Antarctic ice sheet parameters, informed by a previous model calibration to
#              paleo data. There are two options for the Antarctic priors (1) fitting marginal distributions using a
#              kernel density estimation or (2) fitting a multivariate normal distribution that accounts for correlations
#              that emerge during the paleo calibration (note, many of the marginal paleo pdfs are not normally distributed).
#
# Function Arguments:
#
#     joint_antarctic_prior = TRUE/FALSE check for whether to use a joint normal prior distribution (TRUE = option 1 described
#                             above) or fitted marginal kernel density estimates (FLASE = option 2 described above).
#     uniform_ECS           = TRUE/FALSE check for whether or not to use a uniform prior distribution for the equilibrium
#                             climate sensitivity (true = use uniform).
#----------------------------------------------------------------------------------------------------------------------

function construct_sneasybrick_log_prior(joint_antarctic_prior::Bool, uniform_ECS::Bool)

    #---------------------------------------------
    # Antarctic ice sheet priors
    #---------------------------------------------

    # Calculate upper and lower bounds for Antarctic ice sheet parameters (min/max values from paleo calibration).
    antarctic_lower_bound = vec(minimum(antarctic_paleo_params, dims=1))
    antarctic_upper_bound = vec(maximum(antarctic_paleo_params, dims=1))

    # Initialize a vector to store sampled Antarctic ice sheet parameters.
    antarctic_params = zeros(15)

    # Create either the joint normal prior distribution for the Antarctic ice sheet model or the marginal kernel density estimates.
    if joint_antarctic_prior == true

        # Fit a multivariate normal to the data (do not use variance estimate for paleo data, separately estimate AR(1) parameters for recent observations isntead).
        antarctic_joint_prior = fit(MvNormal, antarctic_paleo_params')

    else
        # If not using joint prior, create marginal kernel density estimates for each Antarctic ice sheet model parameter.
        prior_anto_α         = truncated_kernel(antarctic_paleo_params[:,1],  antarctic_lower_bound[1],  antarctic_upper_bound[1])
        prior_anto_β         = truncated_kernel(antarctic_paleo_params[:,2],  antarctic_lower_bound[2],  antarctic_upper_bound[2])
        prior_γ              = truncated_kernel(antarctic_paleo_params[:,3],  antarctic_lower_bound[3],  antarctic_upper_bound[3])
        prior_α              = truncated_kernel(antarctic_paleo_params[:,4],  antarctic_lower_bound[4],  antarctic_upper_bound[4])
        prior_μ              = truncated_kernel(antarctic_paleo_params[:,5],  antarctic_lower_bound[5],  antarctic_upper_bound[5])
        prior_ν              = truncated_kernel(antarctic_paleo_params[:,6],  antarctic_lower_bound[6],  antarctic_upper_bound[6])
        prior_precip₀        = truncated_kernel(antarctic_paleo_params[:,7],  antarctic_lower_bound[7],  antarctic_upper_bound[7])
        prior_κ              = truncated_kernel(antarctic_paleo_params[:,8],  antarctic_lower_bound[8],  antarctic_upper_bound[8])
        prior_flow₀          = truncated_kernel(antarctic_paleo_params[:,9],  antarctic_lower_bound[9],  antarctic_upper_bound[9])
        prior_runoff_height₀ = truncated_kernel(antarctic_paleo_params[:,10], antarctic_lower_bound[10], antarctic_upper_bound[10])
        prior_c              = truncated_kernel(antarctic_paleo_params[:,11], antarctic_lower_bound[11], antarctic_upper_bound[11])
        prior_bedheight₀     = truncated_kernel(antarctic_paleo_params[:,12], antarctic_lower_bound[12], antarctic_upper_bound[12])
        prior_slope          = truncated_kernel(antarctic_paleo_params[:,13], antarctic_lower_bound[13], antarctic_upper_bound[13])
        prior_λ              = truncated_kernel(antarctic_paleo_params[:,14], antarctic_lower_bound[14], antarctic_upper_bound[14])
        prior_temp_threshold = truncated_kernel(antarctic_paleo_params[:,15], antarctic_lower_bound[15], antarctic_upper_bound[15])
    end

    # Create a function to calculate the total log-prior for the Antarctic ice sheet, depending on Antarctic prior type specification.
    antarctic_total_prior =

        # If using the joint normal prior.
        if joint_antarctic_prior == true
            function(p::Array{Float64,1})
                # If sampled parameters are outside paleo range, return -Inf (needed because Distributions.jl package doesn't have truncated multivariate pdfs).
                if any(p .< antarctic_lower_bound) || any(p .> antarctic_upper_bound)
                    return -Inf
                else
                    return logpdf(antarctic_joint_prior, p)
                end
            end

        # Create function using the marginal kernel densities.
        else
            function(p::Array{Float64,1})
                return log(pdf(prior_anto_α, p[1])) + log(pdf(prior_anto_β, p[2])) + log(pdf(prior_γ, p[3])) + log(pdf(prior_α, p[4])) + log(pdf(prior_μ, p[5])) +
                       log(pdf(prior_ν, p[6])) + log(pdf(prior_precip₀, p[7])) + log(pdf(prior_κ, p[8])) + log(pdf(prior_flow₀, p[9])) +
                       log(pdf(prior_runoff_height₀, p[10])) + log(pdf(prior_c, p[11])) + log(pdf(prior_bedheight₀, p[12])) + log(pdf(prior_slope, p[13])) +
                       log(pdf(prior_λ, p[14])) + log(pdf(prior_temp_threshold, p[15]))
            end
        end

    # Declare prior distributions for all other uncertain model and statistical process parameters.

    # -----------------------------------------
    # Statistical Process Priors.
    # -----------------------------------------
    prior_σ_temperature      = Uniform(1e-10, 0.2) # SC-CH4 and Klaus papers.
    prior_σ_ocean_heat       = Uniform(1e-10, 4) # From SC-CH4 paper
    prior_σ_glaciers         = Uniform(1e-10, 0.0015) # Based on BRICK code.
    prior_σ_greenland        = Uniform(1e-10, 0.002) # Based on BRICK code.
    prior_σ_antarctic        = Uniform(1e-10, 0.063) # Based on BRICK code.
    prior_σ_gmsl             = Uniform(1e-10, 0.05) # Just setting the same as prior_σ_gmsl_1900 value from old BRICK code.
    prior_σ²_white_noise_CO₂ = Uniform(1e-10, 200) # Sc-CH4 paper.

    prior_ρ_temperature      = Uniform(-0.99, 0.99)
    prior_ρ_ocean_heat       = Uniform(-0.99, 0.99)
    prior_ρ_glaciers         = Uniform(-0.99, 0.99)
    prior_ρ_greenland        = Uniform(-0.99, 0.99)
    prior_ρ_antarctic        = Uniform(-0.99, 0.99)
    prior_ρ_gmsl             = Truncated(Normal(0.8, .25), -1.0, 1.0)
    prior_α₀_CO₂             = Uniform(0.01, 11.5) # SC-CH4 paper.

    # -----------------------------------------
    # Initial Condition Priors.
    # -----------------------------------------
    prior_CO₂_0              = Uniform(275, 281)
    prior_N₂O_0              = Uniform(264, 282)
    prior_temperature_0      = Truncated(Normal(),-0.3, 0.3)
    prior_ocean_heat_0       = Uniform(-100, 0)
    prior_thermal_s₀         = Uniform(-0.0484, 0.0484) # BRICK defaults. # Initial sea level rise due to thermal expansion designated in 1850 (m SLE).
    prior_greenland_v₀       = Uniform(7.16, 7.56)
    prior_glaciers_v₀        = Uniform(0.31, 0.53)
    prior_glaciers_s₀        = Uniform(-0.0536, 0.0791)
    prior_antarctic_s₀       = Uniform(-0.04755, 0.05585) # Informed by prior BRICK runs.

    # -----------------------------------------
    # Carbon Cycle Priors.
    # -----------------------------------------
    prior_Q10                = Uniform(1.0, 5)
    prior_CO₂_fertilization  = Uniform(0., 1)
    prior_CO₂_diffusivity    = Uniform(0., 200)

    # -----------------------------------------
    # Climate & Radiative Forcing Priors.
    # -----------------------------------------
    prior_heat_diffusivity   = LogNormal(1.1, 0.3)
    #prior_heat_diffusivity   = Uniform(0.1, 4) # from BRICK paper
    prior_rf_scale_aerosol   = TriangularDist(0., 3., 1.)

    # Decide whether to use a uniform or paleo-informed ECS prior.
    if uniform_ECS == true
        prior_ECS = Uniform(0.0, 10.0)
    else
        prior_ECS = Truncated(Cauchy(3.0,2.0), 0.0, 10.0)
    end

    # ---------------------------------------------
    # Sea Level Rise From Thermal Expansion Priors
    # ---------------------------------------------
    prior_thermal_α          = Uniform(0.05, 0.3) # upper/lower bounds from "Impacts of Observational Constraints Related to Sea Level on Estimates of Climate Sensitivity"  # Global ocean-averaged thermal expansion coefficient (kg m⁻³ °C⁻¹).

    #------------------------------------------------
    # Sea Level Rise from Greenland Ice Sheet Priors
    #------------------------------------------------
    prior_greenland_a        = Uniform(-4.0, -0.001)
    prior_greenland_b        = Uniform(5.888, 8.832)
    prior_greenland_α        = Uniform(0.0, 0.001)
    prior_greenland_β        = Uniform(0.0, 0.001)

    #------------------------------------------------------
    # Sea Level Rise from Glaciers & Small Ice Caps Priors
    #------------------------------------------------------
    prior_glaciers_β₀        = Uniform(0.0, 0.041)
    prior_glaciers_n         = Uniform(0.55, 1.0)

    #------------------------------------------------------------------------------------
    # Create function that returns the log-prior of all uncertain model parameters.
    #------------------------------------------------------------------------------------

    function total_log_prior(p::Array{Float64,1})

        # Assign parameter values names for convenience/tractability.
        σ_temperature            = p[1]
        σ_ocean_heat             = p[2]
        σ_glaciers               = p[3]
        σ_greenland              = p[4]
        σ_antarctic              = p[5]
        σ_gmsl                   = p[6]
        σ²_white_noise_CO₂       = p[7]
        ρ_temperature            = p[8]
        ρ_ocean_heat             = p[9]
        ρ_glaciers               = p[10]
        ρ_greenland              = p[11]
        ρ_antarctic              = p[12]
        ρ_gmsl                   = p[13]
        α₀_CO₂                   = p[14]
        CO₂_0                    = p[15]
        N₂O_0                    = p[16]
        temperature_0            = p[17]
        ocean_heat_0             = p[18]
        thermal_s₀               = p[19]
        greenland_v₀             = p[20]
        glaciers_v₀              = p[21]
        glaciers_s₀              = p[22]
        antarctic_s₀             = p[23]
        Q10                      = p[24]
        CO₂_fertilization        = p[25]
        CO₂_diffusivity          = p[26]
        heat_diffusivity         = p[27]
        rf_scale_aerosol         = p[28]
        ECS                      = p[29]
        thermal_α                = p[30]
        greenland_a              = p[31]
        greenland_b              = p[32]
        greenland_α              = p[33]
        greenland_β              = p[34]
        glaciers_β₀              = p[35]
        glaciers_n               = p[36]
        antarctic_params[:]      = p[37:51]

        log_prior = logpdf(prior_σ_temperature, σ_temperature) + logpdf(prior_σ_ocean_heat, σ_ocean_heat) + logpdf(prior_σ_glaciers, σ_glaciers) + logpdf(prior_σ_greenland, σ_greenland) + logpdf(prior_σ_antarctic, σ_antarctic) + logpdf(prior_σ_gmsl, σ_gmsl) + logpdf(prior_σ²_white_noise_CO₂, σ²_white_noise_CO₂) +
                    logpdf(prior_ρ_temperature, ρ_temperature) + logpdf(prior_ρ_ocean_heat, ρ_ocean_heat) + logpdf(prior_ρ_glaciers, ρ_glaciers) + logpdf(prior_ρ_greenland, ρ_greenland) + logpdf(prior_ρ_antarctic, ρ_antarctic) + logpdf(prior_ρ_gmsl, ρ_gmsl) + logpdf(prior_α₀_CO₂, α₀_CO₂) +
                    logpdf(prior_CO₂_0, CO₂_0) + logpdf(prior_N₂O_0, N₂O_0) + logpdf(prior_temperature_0, temperature_0) + logpdf(prior_ocean_heat_0, ocean_heat_0) + logpdf(prior_thermal_s₀, thermal_s₀) + logpdf(prior_greenland_v₀, greenland_v₀) + logpdf(prior_glaciers_v₀, glaciers_v₀) + logpdf(prior_glaciers_s₀, glaciers_s₀) + logpdf(prior_antarctic_s₀, antarctic_s₀) +
                    logpdf(prior_Q10, Q10) + logpdf(prior_CO₂_fertilization, CO₂_fertilization) + logpdf(prior_CO₂_diffusivity, CO₂_diffusivity) +
                    logpdf(prior_heat_diffusivity, heat_diffusivity) + logpdf(prior_rf_scale_aerosol, rf_scale_aerosol) + logpdf(prior_ECS, ECS) +
                    logpdf(prior_thermal_α, thermal_α) +
                    logpdf(prior_greenland_a, greenland_a) + logpdf(prior_greenland_b, greenland_b) + logpdf(prior_greenland_α, greenland_α) + logpdf(prior_greenland_β, greenland_β) +
                    logpdf(prior_glaciers_β₀, glaciers_β₀) + logpdf(prior_glaciers_n, glaciers_n) +
                    antarctic_total_prior(antarctic_params)

        return log_prior
    end

    # Return total log-prior function for SNEASY+BRICK.
    return total_log_prior
end



#######################################################################################################################
# CALCULATE LOG POSTERIOR.
#######################################################################################################################
# Description: This creates a function that calculates the log-posterior probability of the uncertain model, initial
#              condition, and statistical process parameters.
#
# Function Arguments:
#
#     f_run_model           = A function that runs the specific climate model version and returns the output being calibrated to observations.
#     model_start_year      = First year to run the model (not necessarily first year of the calibration if model initializes earlier).
#     end_year              = The final year to run the model calibration (defaults to 2017).
#     joint_antarctic_prior = TRUE/FALSE check for whether to use a joint normal prior distribution (TRUE = option 1 described
#                             above) or fitted marginal kernel density estimates (FLASE = option 2 described above).
#     uniform_ECS           = TRUE/FALSE check for whether or not to use a uniform prior distribution for the equilibrium
#                             climate sensitivity (true = use uniform).
#----------------------------------------------------------------------------------------------------------------------

function construct_sneasybrick_log_posterior(f_run_model!; model_start_year::Int=1850, calibration_end_year::Int=2017, joint_antarctic_prior::Bool=false, uniform_ECS::Bool=false)

   # Create a vector of calibration years and calculate total number of years to run model.
    calibration_years = collect(model_start_year:calibration_end_year)
    n = length(calibration_years)

    # Get log-prior function.
    sneasybrick_log_prior = construct_sneasybrick_log_prior(joint_antarctic_prior, uniform_ECS)

    # Load calibration data/observations.
    calibration_data, obs_antarctic_trends, obs_thermal_trends = load_calibration_data(model_start_year, calibration_end_year, last_sea_level_norm_year=1990)

    # Calculate indices for each year that has an observation in calibration data sets.
    indices_maunaloa_co2_data  = findall(x-> !ismissing(x), calibration_data.maunaloa_co2_obs)
    indices_lawdome_co2_data   = findall(x-> !ismissing(x), calibration_data.lawdome_co2_obs)
    indices_oceanco2_flux_data = findall(x-> !ismissing(x), calibration_data.oceanco2_flux_obs)
    indices_temperature_data   = findall(x-> !ismissing(x), calibration_data.hadcrut_temperature_obs)
    indices_oceanheat_data     = findall(x-> !ismissing(x), calibration_data.ocean_heat_obs)
    indices_glaciers_data      = findall(x-> !ismissing(x), calibration_data.glaciers_obs)
    indices_greenland_data     = findall(x-> !ismissing(x), calibration_data.merged_greenland_obs) # Use merged Greenland data.
    indices_antarctic_data     = findall(x-> !ismissing(x), calibration_data.antarctic_imbie_obs)
    indices_gmsl_data          = findall(x-> !ismissing(x), calibration_data.gmsl_obs)

    # Combine CO₂ indices from Law Dome and Mauna Loa observations.
    indices_co2_data = sort(vcat(indices_lawdome_co2_data, indices_maunaloa_co2_data))

    # Combine CO₂ measurement errors from Law Dome and Mauna Loa observations (just for convenience).
    calibration_data.co2_combined_sigma = calibration_data.lawdome_co2_sigma
    calibration_data.co2_combined_sigma[indices_maunaloa_co2_data] = calibration_data.maunaloa_co2_sigma[indices_maunaloa_co2_data]

    # Calculate number of ice core observations for CO₂ (used for indexing).
    n_lawdome_co2 = length(indices_lawdome_co2_data)

    # Allocate arrays to store data-model residuals.
    co2_residual         = zeros(length(indices_co2_data))
    temperature_residual = zeros(length(indices_temperature_data))
    ocean_heat_residual  = zeros(length(indices_oceanheat_data))
    glaciers_residual    = zeros(length(indices_glaciers_data))
    greenland_residual   = zeros(length(indices_greenland_data))
    antarctic_residual   = zeros(length(indices_antarctic_data))
    gmsl_residual        = zeros(length(indices_gmsl_data))

    # Allocate vectors to store model output being calibrated to the observations.
    modeled_CO₂               = zeros(n)
    modeled_oceanCO₂_flux     = zeros(n)
    modeled_temperature       = zeros(n)
    modeled_ocean_heat        = zeros(n)
    modeled_glaciers          = zeros(n)
    modeled_greenland         = zeros(n)
    modeled_antarctic         = zeros(n)
    modeled_thermal_expansion = zeros(n)
    modeled_thermal_trend     = zeros(size(obs_thermal_trends, 1))
    modeled_gmsl              = zeros(n)

    # Allocate vectors to store log-likelihoods for individual thermal trends and ocean CO₂ flux data points for convenience (assumes iid error structure).
    individual_llik_thermal_trend = zeros(size(obs_thermal_trends, 1))
    individual_llik_oceanco2_flux = zeros(length(indices_oceanco2_flux_data))

    #---------------------------------------------------------------------------------------------------------------------------------------
    # Create a function to calculate the log-likelihood for the observations, assuming residual independence across calibration data sets.
    #---------------------------------------------------------------------------------------------------------------------------------------

    function sneasybrick_log_likelihood(p::Array{Float64,1})

        # Assign names to uncertain statistical process parameters used in log-likelihood calculations.
        σ_temperature      = p[1]
        σ_ocean_heat       = p[2]
        σ_glaciers         = p[3]
        σ_greenland        = p[4]
        σ_antarctic        = p[5]
        σ_gmsl             = p[6]
        σ²_white_noise_CO₂ = p[7]
        ρ_temperature      = p[8]
        ρ_ocean_heat       = p[9]
        ρ_glaciers         = p[10]
        ρ_greenland        = p[11]
        ρ_antarctic        = p[12]
        ρ_gmsl             = p[13]
        α₀_CO₂             = p[14]

        # Run an instance of SNEASY+BRICK with sampled parameter set and return model output being compared to observations.
        f_run_model!(p, modeled_CO₂, modeled_oceanCO₂_flux, modeled_temperature, modeled_ocean_heat,
                     modeled_glaciers, modeled_greenland, modeled_antarctic, modeled_thermal_expansion, modeled_gmsl)

        #---------------------------------------------------------------------------
        # Global Surface Temperature (normalized to 1861-1880 mean) Log-Likelihood.
        #---------------------------------------------------------------------------

        llik_temperature = 0.0

        # Calculate temperature residuals.
        for (i, index)=enumerate(indices_temperature_data)
            temperature_residual[i] = calibration_data[index, :hadcrut_temperature_obs] - modeled_temperature[index]
        end

        # Calculate temperature log-likelihood.
        llik_temperature = hetero_logl_ar1(temperature_residual, σ_temperature, ρ_temperature, calibration_data[indices_temperature_data,:hadcrut_temperature_sigma])

        #-----------------------------------------------------------------------
        # Ocean Heat Content Log-Likelihood
        #-----------------------------------------------------------------------

        llik_ocean_heat = 0.0

        # Calculate ocean heat residuals.
        for (i, index)=enumerate(indices_oceanheat_data)
            ocean_heat_residual[i] = calibration_data[index, :ocean_heat_obs] - modeled_ocean_heat[index]
        end

        # Calculate ocean heat log-likelihood.
        llik_ocean_heat = hetero_logl_ar1(ocean_heat_residual, σ_ocean_heat, ρ_ocean_heat, calibration_data[indices_oceanheat_data, :ocean_heat_sigma])

        #-----------------------------------------------------------------------
        # Atmospheric CO₂ Concentration Log-Likelihood
        #-----------------------------------------------------------------------

        llik_co2 = 0.

        # Calculate CO₂ concentration (Law Dome) residuals (assuming 8 year model mean centered on year of ice core observation).
        for (i, index)=enumerate(indices_lawdome_co2_data)
            co2_residual[i] = calibration_data[index, :lawdome_co2_obs] - mean(modeled_CO₂[index .+ (-4:3)])
        end

        # Calculate CO₂ concentration (Mauna Loa) residuals.
        for (i, index)=enumerate(indices_maunaloa_co2_data)
            co2_residual[i+n_lawdome_co2] = calibration_data[index, :maunaloa_co2_obs] - modeled_CO₂[index]
        end

        # Calculate atmospheric CO₂ concentration log-likelihood.
        llik_co2 = hetero_logl_car1(co2_residual, indices_co2_data, σ²_white_noise_CO₂, α₀_CO₂, calibration_data[indices_co2_data, :co2_combined_sigma])

        #-----------------------------------------------------------------------
        # Ocean Carbon Flux Log-Likelihood
        #-----------------------------------------------------------------------

        llik_oceanco2_flux = 0.

        # Calculate ocean CO₂ flux log-likelihood for individual data points.
        for (i,index) = enumerate(indices_oceanco2_flux_data)
            individual_llik_oceanco2_flux[i] = logpdf(Normal(modeled_oceanCO₂_flux[index], calibration_data[index, :oceanco2_flux_sigma]), calibration_data[index, :oceanco2_flux_obs])
        end

        # Calculate ocean CO₂ flux total log-likelihood as sum of individual data point likelihoods.
        llik_oceanco2_flux = sum(individual_llik_oceanco2_flux)

        #-----------------------------------------------------------------------
        # Glaciers and Small Ice Caps Log-Likelihood
        #-----------------------------------------------------------------------

        llik_glaciers = 0.0

        # Calculate glaciers and small ice caps residuals.
        for (i, index)=enumerate(indices_glaciers_data)
            glaciers_residual[i] = calibration_data[index, :glaciers_obs] - modeled_glaciers[index]
        end

        # Calculate sea level contribution from glaciers and small ice caps log-likelihood.
        llik_glaciers = hetero_logl_ar1(glaciers_residual, σ_glaciers, ρ_glaciers, calibration_data[indices_glaciers_data, :glaciers_sigma])

        #-------------------------------------------------------------------------------
        # Greenland Ice Sheet Merged Data (normalized to 1992_2001 mean) Log-Likelihood
        #-------------------------------------------------------------------------------

        llik_greenland = 0.0

        # Calculate greenland ice sheet residuals.
        for (i, index)=enumerate(indices_greenland_data)
            greenland_residual[i] = calibration_data[index, :merged_greenland_obs] - modeled_greenland[index]
        end

        # Calculate sea level contribution from Greenland ice sheet log-likelihood.
        llik_greenland = hetero_logl_ar1(greenland_residual, σ_greenland, ρ_greenland, calibration_data[indices_greenland_data, :merged_greenland_sigma])

        #------------------------------------------------------------------------------------
        # AIS (Antarctic Ice Sheet) IMBIE Data (normalized to 1992_2001 mean) Log-Likelihood
        #------------------------------------------------------------------------------------

        llik_antarctic = 0.0

        # Calculate Antarctic ice sheet residuals.
        for (i, index)=enumerate(indices_antarctic_data)
            antarctic_residual[i] = calibration_data[index, :antarctic_imbie_obs] - modeled_antarctic[index]
        end

        # Calculate sea level contribution from Antarctic ice sheet log-likelihood.
        llik_antarctic = hetero_logl_ar1(antarctic_residual, σ_antarctic, ρ_antarctic, calibration_data[indices_antarctic_data, :antarctic_imbie_sigma])

        #-----------------------------------------------------------------------
        # Thermal Expansion Trends for 1971-2009 & 1993-2009
        #-----------------------------------------------------------------------
        llik_thermal_trend = 0.0

        # Add a check for whether or not the calibration period covers the thermal trend period (some sensitivity tests will not), otherwise just return 0.0.
        if calibration_end_year >= 2009

            # Calculate the AIS trends (in milimeters) from the annual modeled output.
            modeled_thermal_trend[:] = calculate_trends(modeled_thermal_expansion, obs_thermal_trends, model_start_year, calibration_end_year)

            # Calculate thermal expansion trend residuals.
            for i = 1:length(modeled_thermal_trend)
                thermal_trend_residual = modeled_thermal_trend[i] - obs_thermal_trends.Trend[i]
                # Calculate σ for observed trends based on IPCC 90% trend window values.
                obs_trend_err = 0.5 * (obs_thermal_trends.Upper_90_Percent[i] - obs_thermal_trends.Lower_90_Percent[i])
                # Calculate sea level contribution from thermal expansion log-likelihood.
                individual_llik_thermal_trend[i] = logpdf(Normal(0.0, obs_trend_err), thermal_trend_residual)
            end

            # Calculate total log-likelihood of thermal expansion sea level trend as sum of individual log likelihoods.
            llik_thermal_trend = sum(individual_llik_thermal_trend)
        end

        #---------------------------------------------------------------------------
        # Global Mean Sea Level Rise (normalized to 1961-1990 mean) Log-Likelihood.
        #---------------------------------------------------------------------------

        llik_gmsl = 0.0

        # Calculate global mean sea level residuals.
        for (i, index)=enumerate(indices_gmsl_data)
            gmsl_residual[i] = calibration_data[index, :gmsl_obs] - modeled_gmsl[index]
        end

        # Calculate global mean sea level log-likelihood.
        llik_gmsl = hetero_logl_ar1(gmsl_residual, σ_gmsl, ρ_gmsl, calibration_data[indices_gmsl_data,:gmsl_sigma])

        #-----------------------------------------------------------------------
        # Total Log-Likelihood
        #-----------------------------------------------------------------------

        # Calculate the total log-likelihood (assuming residual independence across data sets).
        llik = llik_temperature + llik_ocean_heat + llik_co2 + llik_oceanco2_flux + llik_glaciers + llik_greenland + llik_antarctic + llik_thermal_trend + llik_gmsl

        return llik
    end

    #---------------------------------------------------------------------------------------------------------------
    # Create a function to calculate the log-posterior of uncertain parameters 'p' (posterior ∝ likelihood * prior)
    #---------------------------------------------------------------------------------------------------------------

    function sneasybrick_log_posterior(p)

        # Calculate log-prior
        log_prior = sneasybrick_log_prior(p)

        # In case a parameter sample leads to non-physical model outcomes, return -Inf rather than erroring out.
        try
            log_post = isfinite(log_prior) ? sneasybrick_log_likelihood(p) + log_prior : -Inf
        catch
            log_post = - Inf
        end
    end

    # Return log posterior function given user specifications.
    return sneasybrick_log_posterior
end
