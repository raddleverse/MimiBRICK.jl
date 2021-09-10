using CSVFiles
using DataFrames
using Distributions
using Mimi
using Statistics

# Load files with additional functions needed for analysis.
include(joinpath("..", "helper_functions.jl"))
include(joinpath("..", "..", "calibration", "calibration_helper_functions.jl"))

# Load the default version of SNEASY+BRICK.
include(joinpath("..", "create_models", "SNEASY_BRICK.jl"))

#############################################################################################################################
# RUN BASELINE VERSION OF SNEASY+BRICK.
#############################################################################################################################
# Description: This creates a function that runs two baseline versions of SNEASY+BRICK (a standard run, and a run with
#              an extra pulse of CO₂ emissions in a user-specified year) and saves the key model projection output.
#
# Function Arguments:
#
#       rcp_scenario  = A string identifying the RCP emissions and forcing scenario to use (options are "RCP26" and "RCP85").
#       pulse_year    = The year to add a pulse of methane emissions.
#       pulse_size    = The size of the methane emissions pulse in MtCH₄.
#       end_year      = The final year to run the model for.
#----------------------------------------------------------------------------------------------------------------------------

function construct_sneasybrick_baseline(rcp_scenario::String,  pulse_year::Int, pulse_size::Float64, end_year::Int)

    #------------------------------------------
    # Set up emissions scenario and CO₂ pulse.
    #------------------------------------------

    # Set first year to run model.
    start_year = 1850

    # Set model years and calculate number of years.
    model_years  = collect(start_year:end_year)
    number_years = length(model_years)

    # Find indices for RCP data (1765-2500) corresponding to DICE years.
    rcp_indices = findall((in)(model_years), 1765:2500)

    # Load emissions data (index into appropriate years).
    rcp_emissions = DataFrame(load(joinpath(@__DIR__, "..", "..", "data", "model_data", rcp_scenario*"_emissions.csv"), skiplines_begin=36))

    # Set up baseline CO₂ emissions time series that will have the pulse.
    rcp_co2_emissions_pulse = (rcp_emissions.FossilCO2 .+ rcp_emissions.OtherCO2)[rcp_indices]

    # Add the emission pulse in user-specified year.
    pulse_year_index = findfirst(x -> x == pulse_year, model_years)
    rcp_co2_emissions_pulse[pulse_year_index] = rcp_co2_emissions_pulse[pulse_year_index] + pulse_size

    #------------------------------------------
    # Set up SNEASY+BRICK runs.
    #------------------------------------------

     # Get indices needed to normalize temperature anomalies relative to 1861-1880 mean.
    temperature_norm_indices = findall((in)(1861:1880), start_year:end_year)

    # Get indices needed to normalize all sea level rise sources.
    sealevel_norm_indices_1961_1990 = findall((in)(1961:1990), start_year:end_year)
    sealevel_norm_indices_1992_2001 = findall((in)(1992:2001), start_year:end_year)

    # Load calibration data from 1765-2017 (measurement errors used in simulated noise).
    calibration_data, obs_antarctic_trends, obs_thermal_trends = load_calibration_data(start_year, 2017)

    # Pre-allocate vectors to hold simulated CAR(1) & AR(1) with measurement error noise.
    normal_noise_oceanco2 = zeros(number_years)
    car1_noise_co2        = zeros(number_years)
    ar1_noise_temperature = zeros(number_years)
    ar1_noise_ocean_heat  = zeros(number_years)
    ar1_noise_glaciers    = zeros(number_years)
    ar1_noise_greenland   = zeros(number_years)
    ar1_noise_antarctic   = zeros(number_years)
    ar1_noise_gmsl        = zeros(number_years)

    # Replicate errors for years without observations over model time horizon (used for simulating AR1 noise).
    obs_error_temperature = replicate_errors(start_year, end_year, calibration_data.hadcrut_temperature_sigma)
    obs_error_oceanheat   = replicate_errors(start_year, end_year, calibration_data.ocean_heat_sigma)
    obs_error_glaciers    = replicate_errors(start_year, end_year, calibration_data.glaciers_sigma)
    obs_error_greenland   = replicate_errors(start_year, end_year, calibration_data.merged_greenland_sigma)
    obs_error_antarctic   = replicate_errors(start_year, end_year, calibration_data.antarctic_imbie_sigma)
    obs_error_gmsl        = replicate_errors(start_year, end_year, calibration_data.gmsl_sigma)

    # For CAR(1) noise, set constant CO₂ observations errors for start year-1958 (Law Dome), and 1959-end year (Mauna Loa).
    obs_error_co2 = ones(number_years) .* unique(skipmissing(calibration_data.maunaloa_co2_sigma))[1]
    obs_error_co2[1:findfirst(isequal(1958), model_years)] .= unique(skipmissing(calibration_data.lawdome_co2_sigma))[1]

    # Also need to calculate landwater storage contribution so it is the same between base and pulse runs.
    landwater_storage_sl = zeros(number_years)

    # Create two versions of SNEASY+BRICK (base model and pulse model)
    sneasybrick_base  = create_sneasy_brick(rcp_scenario, end_year = end_year)
    sneasybrick_pulse = create_sneasy_brick(rcp_scenario, end_year = end_year)

    # Update pulse model to use emission scneario with an extra emissions pulse in a single year.
    update_param!(sneasybrick_pulse, :ccm, :CO2_emissions, rcp_co2_emissions_pulse)

    #---------------------------------------------------------------------------------------------------------------------
    # Create a function to run SNEASY+BRICK over the calibrated posterior model parameters.
    #---------------------------------------------------------------------------------------------------------------------

    function sneasybrick_baseline(calibrated_parameters::Array{Float64,2}, ci_interval_1::Float64, ci_interval_2::Float64)

        # Calculate number of calibrated parameter samples (row = sample from joint posterior distribution, column = specific parameter).
        number_samples = size(calibrated_parameters, 1)

        # Pre-allocate arrays to store SNEASY+BRICK results.
        base_temperature  = zeros(Union{Missing, Float64}, number_samples, number_years)
        base_co2          = zeros(Union{Missing, Float64}, number_samples, number_years)
        base_ocean_heat   = zeros(Union{Missing, Float64}, number_samples, number_years)
        base_oceanco2     = zeros(Union{Missing, Float64}, number_samples, number_years)
        base_glaciers     = zeros(Union{Missing, Float64}, number_samples, number_years)
        base_greenland    = zeros(Union{Missing, Float64}, number_samples, number_years)
        base_thermal_sl   = zeros(Union{Missing, Float64}, number_samples, number_years)
        base_antarctic    = zeros(Union{Missing, Float64}, number_samples, number_years)
        base_gmsl         = zeros(Union{Missing, Float64}, number_samples, number_years)

        pulse_temperature = zeros(Union{Missing, Float64}, number_samples, number_years)
        pulse_co2         = zeros(Union{Missing, Float64}, number_samples, number_years)
        pulse_gmsl        = zeros(Union{Missing, Float64}, number_samples, number_years)

        # For each calibrated parameter sample, run the base and pulse versions of SNEASY+BRICK, superimpose noise, and store results.
        for i in 1:number_samples

            # Assign sampled parameters names to keep things organized.
            σ_temperature            = calibrated_parameters[i,1]
            σ_ocean_heat             = calibrated_parameters[i,2]
            σ_glaciers               = calibrated_parameters[i,3]
            σ_greenland              = calibrated_parameters[i,4]
            σ_antarctic              = calibrated_parameters[i,5]
            σ_gmsl                   = calibrated_parameters[i,6]
            σ²_white_noise_CO₂       = calibrated_parameters[i,7]
            ρ_temperature            = calibrated_parameters[i,8]
            ρ_ocean_heat             = calibrated_parameters[i,9]
            ρ_glaciers               = calibrated_parameters[i,10]
            ρ_greenland              = calibrated_parameters[i,11]
            ρ_antarctic              = calibrated_parameters[i,12]
            ρ_gmsl                   = calibrated_parameters[i,13]
            α₀_CO₂                   = calibrated_parameters[i,14]
            CO₂_0                    = calibrated_parameters[i,15]
            N₂O_0                    = calibrated_parameters[i,16]
            temperature_0            = calibrated_parameters[i,17]
            ocean_heat_0             = calibrated_parameters[i,18]
            thermal_s₀               = calibrated_parameters[i,19]
            greenland_v₀             = calibrated_parameters[i,20]
            glaciers_v₀              = calibrated_parameters[i,21]
            glaciers_s₀              = calibrated_parameters[i,22]
            antarctic_s₀             = calibrated_parameters[i,23]
            Q10                      = calibrated_parameters[i,24]
            CO₂_fertilization        = calibrated_parameters[i,25]
            CO₂_diffusivity          = calibrated_parameters[i,26]
            heat_diffusivity         = calibrated_parameters[i,27]
            rf_scale_aerosol         = calibrated_parameters[i,28]
            ECS                      = calibrated_parameters[i,29]
            thermal_α                = calibrated_parameters[i,30]
            greenland_a              = calibrated_parameters[i,31]
            greenland_b              = calibrated_parameters[i,32]
            greenland_α              = calibrated_parameters[i,33]
            greenland_β              = calibrated_parameters[i,34]
            glaciers_β₀              = calibrated_parameters[i,35]
            glaciers_n               = calibrated_parameters[i,36]
            anto_α                   = calibrated_parameters[i,37]
            anto_β                   = calibrated_parameters[i,38]
            antarctic_γ              = calibrated_parameters[i,39]
            antarctic_α              = calibrated_parameters[i,40]
            antarctic_μ              = calibrated_parameters[i,41]
            antarctic_ν              = calibrated_parameters[i,42]
            antarctic_precip₀        = calibrated_parameters[i,43]
            antarctic_κ              = calibrated_parameters[i,44]
            antarctic_flow₀          = calibrated_parameters[i,45]
            antarctic_runoff_height₀ = calibrated_parameters[i,46]
            antarctic_c              = calibrated_parameters[i,47]
            antarctic_bedheight₀     = calibrated_parameters[i,48]
            antarctic_slope          = calibrated_parameters[i,49]
            antarctic_λ              = calibrated_parameters[i,50]
            antarctic_temp_threshold = calibrated_parameters[i,51]

            # Set parameters for base version of SNEASY+BRICK.

            update_param!(sneasybrick_base, :model_CO₂_0, CO₂_0) # shared parameter linked to both :rfco2 and :ccm

            update_param!(sneasybrick_base, :doeclim, :t2co, ECS)
            update_param!(sneasybrick_base, :doeclim, :kappa, heat_diffusivity)

            update_param!(sneasybrick_base, :ccm, :Q10, Q10)
            update_param!(sneasybrick_base, :ccm, :Beta, CO₂_fertilization)
            update_param!(sneasybrick_base, :ccm, :Eta, CO₂_diffusivity)

            update_param!(sneasybrick_base, :rfco2, :N₂O_0, N₂O_0)

            update_param!(sneasybrick_base, :radiativeforcing, :alpha, rf_scale_aerosol)

            update_param!(sneasybrick_base, :antarctic_ocean, :anto_α, anto_α)
            update_param!(sneasybrick_base, :antarctic_ocean, :anto_β, anto_β)

            update_param!(sneasybrick_base, :antarctic_icesheet, :ais_sea_level₀, antarctic_s₀)
            update_param!(sneasybrick_base, :antarctic_icesheet, :ais_bedheight₀, antarctic_bedheight₀)
            update_param!(sneasybrick_base, :antarctic_icesheet, :ais_slope, antarctic_slope)
            update_param!(sneasybrick_base, :antarctic_icesheet, :ais_μ, antarctic_μ)
            update_param!(sneasybrick_base, :antarctic_icesheet, :ais_runoffline_snowheight₀, antarctic_runoff_height₀)
            update_param!(sneasybrick_base, :antarctic_icesheet, :ais_c, antarctic_c)
            update_param!(sneasybrick_base, :antarctic_icesheet, :ais_precipitation₀, antarctic_precip₀)
            update_param!(sneasybrick_base, :antarctic_icesheet, :ais_κ, antarctic_κ)
            update_param!(sneasybrick_base, :antarctic_icesheet, :ais_ν, antarctic_ν)
            update_param!(sneasybrick_base, :antarctic_icesheet, :ais_iceflow₀, antarctic_flow₀)
            update_param!(sneasybrick_base, :antarctic_icesheet, :ais_γ, antarctic_γ)
            update_param!(sneasybrick_base, :antarctic_icesheet, :ais_α, antarctic_α)
            update_param!(sneasybrick_base, :antarctic_icesheet, :temperature_threshold, antarctic_temp_threshold)
            update_param!(sneasybrick_base, :antarctic_icesheet, :λ, antarctic_λ)

            update_param!(sneasybrick_base, :glaciers_small_icecaps, :gsic_β₀, glaciers_β₀)
            update_param!(sneasybrick_base, :glaciers_small_icecaps, :gsic_v₀, glaciers_v₀)
            update_param!(sneasybrick_base, :glaciers_small_icecaps, :gsic_s₀, glaciers_s₀)
            update_param!(sneasybrick_base, :glaciers_small_icecaps, :gsic_n, glaciers_n)

            update_param!(sneasybrick_base, :greenland_icesheet, :greenland_a, greenland_a)
            update_param!(sneasybrick_base, :greenland_icesheet, :greenland_b, greenland_b)
            update_param!(sneasybrick_base, :greenland_icesheet, :greenland_α, greenland_α)
            update_param!(sneasybrick_base, :greenland_icesheet, :greenland_β, greenland_β)
            update_param!(sneasybrick_base, :greenland_icesheet, :greenland_v₀, greenland_v₀)

            update_param!(sneasybrick_base, :thermal_expansion, :te_α, thermal_α)
            update_param!(sneasybrick_base, :thermal_expansion, :te_s₀, thermal_s₀)

            # Set parameters for pulse version of SNEASY+BRICK.
            
            update_param!(sneasybrick_pulse, :model_CO₂_0, CO₂_0) # shared parameter linked to both :rfco2 and :ccm

            update_param!(sneasybrick_pulse, :doeclim, :t2co, ECS)
            update_param!(sneasybrick_pulse, :doeclim, :kappa, heat_diffusivity)

            update_param!(sneasybrick_pulse, :ccm, :Q10, Q10)
            update_param!(sneasybrick_pulse, :ccm, :Beta, CO₂_fertilization)
            update_param!(sneasybrick_pulse, :ccm, :Eta, CO₂_diffusivity)

            update_param!(sneasybrick_pulse, :rfco2, :N₂O_0, N₂O_0)

            update_param!(sneasybrick_pulse, :radiativeforcing, :alpha, rf_scale_aerosol)

            update_param!(sneasybrick_pulse, :antarctic_ocean, :anto_α, anto_α)
            update_param!(sneasybrick_pulse, :antarctic_ocean, :anto_β, anto_β)

            update_param!(sneasybrick_pulse, :antarctic_icesheet, :ais_sea_level₀, antarctic_s₀)
            update_param!(sneasybrick_pulse, :antarctic_icesheet, :ais_bedheight₀, antarctic_bedheight₀)
            update_param!(sneasybrick_pulse, :antarctic_icesheet, :ais_slope, antarctic_slope)
            update_param!(sneasybrick_pulse, :antarctic_icesheet, :ais_μ, antarctic_μ)
            update_param!(sneasybrick_pulse, :antarctic_icesheet, :ais_runoffline_snowheight₀, antarctic_runoff_height₀)
            update_param!(sneasybrick_pulse, :antarctic_icesheet, :ais_c, antarctic_c)
            update_param!(sneasybrick_pulse, :antarctic_icesheet, :ais_precipitation₀, antarctic_precip₀)
            update_param!(sneasybrick_pulse, :antarctic_icesheet, :ais_κ, antarctic_κ)
            update_param!(sneasybrick_pulse, :antarctic_icesheet, :ais_ν, antarctic_ν)
            update_param!(sneasybrick_pulse, :antarctic_icesheet, :ais_iceflow₀, antarctic_flow₀)
            update_param!(sneasybrick_pulse, :antarctic_icesheet, :ais_γ, antarctic_γ)
            update_param!(sneasybrick_pulse, :antarctic_icesheet, :ais_α, antarctic_α)
            update_param!(sneasybrick_pulse, :antarctic_icesheet, :temperature_threshold, antarctic_temp_threshold)
            update_param!(sneasybrick_pulse, :antarctic_icesheet, :λ, antarctic_λ)

            update_param!(sneasybrick_pulse, :glaciers_small_icecaps, :gsic_β₀, glaciers_β₀)
            update_param!(sneasybrick_pulse, :glaciers_small_icecaps, :gsic_v₀, glaciers_v₀)
            update_param!(sneasybrick_pulse, :glaciers_small_icecaps, :gsic_s₀, glaciers_s₀)
            update_param!(sneasybrick_pulse, :glaciers_small_icecaps, :gsic_n, glaciers_n)

            update_param!(sneasybrick_pulse, :greenland_icesheet, :greenland_a, greenland_a)
            update_param!(sneasybrick_pulse, :greenland_icesheet, :greenland_b, greenland_b)
            update_param!(sneasybrick_pulse, :greenland_icesheet, :greenland_α, greenland_α)
            update_param!(sneasybrick_pulse, :greenland_icesheet, :greenland_β, greenland_β)
            update_param!(sneasybrick_pulse, :greenland_icesheet, :greenland_v₀, greenland_v₀)

            update_param!(sneasybrick_pulse, :thermal_expansion, :te_α, thermal_α)
            update_param!(sneasybrick_pulse, :thermal_expansion, :te_s₀, thermal_s₀)

            # Calculate land water storage contribution to sea level rise (sampled from Normal distribution) and set same scenario for base and pulse runs.
            landwater_storage_sl[:] = rand(Normal(0.0003, 0.00018), number_years)
            update_param!(sneasybrick_base,  :landwater_storage, :lws_random_sample, landwater_storage_sl)
            update_param!(sneasybrick_pulse, :landwater_storage, :lws_random_sample, landwater_storage_sl)

            # Add a check for cases where non-physical parameter samples cause a model error (may occur during sensitivity experiements).
            try

                # Run both models.
                run(sneasybrick_base)
                run(sneasybrick_pulse)

                # Create noise to superimpose on results using calibrated statistical parameters and measurement noise (note: both models use same estimated noise for each sample).
                ar1_noise_temperature[:] = simulate_ar1_noise(number_years, σ_temperature, ρ_temperature, obs_error_temperature)
                ar1_noise_ocean_heat[:]  = simulate_ar1_noise(number_years, σ_ocean_heat,  ρ_ocean_heat,  obs_error_oceanheat)
                ar1_noise_glaciers[:]    = simulate_ar1_noise(number_years, σ_glaciers,    ρ_glaciers,    obs_error_glaciers)
                ar1_noise_greenland[:]   = simulate_ar1_noise(number_years, σ_greenland,   ρ_greenland,   obs_error_greenland)
                ar1_noise_antarctic[:]   = simulate_ar1_noise(number_years, σ_antarctic,   ρ_antarctic,   obs_error_antarctic)
                ar1_noise_gmsl[:]        = simulate_ar1_noise(number_years, σ_gmsl,        ρ_gmsl,        obs_error_gmsl)
                normal_noise_oceanco2[:] = rand(Normal(0,0.4*sqrt(10)), number_years)

                # CO₂ uses CAR(1) statistical process parameters.
                car1_noise_co2[:] = simulate_car1_noise(number_years, α₀_CO₂, σ²_white_noise_CO₂, obs_error_co2)

                # Store model projections resulting from parameter sample `i` for base model and normalize relative to appropriate time period.
                base_temperature[i,:] = sneasybrick_base[:doeclim, :temp] .- mean(sneasybrick_base[:doeclim, :temp][temperature_norm_indices]) .+ ar1_noise_temperature .+ temperature_0
                base_co2[i,:]         = sneasybrick_base[:ccm, :atmco2] .+ car1_noise_co2
                base_ocean_heat[i,:]  = sneasybrick_base[:doeclim, :heat_mixed] .+ sneasybrick_base[:doeclim, :heat_interior] .+ ar1_noise_ocean_heat .+ ocean_heat_0
                base_oceanco2[i,:]    = sneasybrick_base[:ccm, :atm_oc_flux] .+ normal_noise_oceanco2
                base_glaciers[i,:]    = sneasybrick_base[:glaciers_small_icecaps, :gsic_sea_level] .- mean(sneasybrick_base[:glaciers_small_icecaps, :gsic_sea_level][sealevel_norm_indices_1961_1990]) .+ ar1_noise_glaciers
                base_greenland[i,:]   = sneasybrick_base[:greenland_icesheet, :greenland_sea_level] .- mean(sneasybrick_base[:greenland_icesheet, :greenland_sea_level][sealevel_norm_indices_1992_2001]) .+ ar1_noise_greenland
                base_antarctic[i,:]   = sneasybrick_base[:antarctic_icesheet, :ais_sea_level] .- mean(sneasybrick_base[:antarctic_icesheet, :ais_sea_level][sealevel_norm_indices_1992_2001]) .+ ar1_noise_antarctic
                base_thermal_sl[i,:]  = sneasybrick_base[:thermal_expansion, :te_sea_level]
                base_gmsl[i,:]        = sneasybrick_base[:global_sea_level, :sea_level_rise] .- mean(sneasybrick_base[:global_sea_level, :sea_level_rise][sealevel_norm_indices_1961_1990]) .+ ar1_noise_gmsl

                # Store and normalize tempeature, CO₂, and sea level projections resulting from parameter sample `i` for pulse model (used for estimating the SC-CO₂).
                pulse_temperature[i,:] = sneasybrick_pulse[:doeclim, :temp] .- mean(sneasybrick_pulse[:doeclim, :temp][temperature_norm_indices]) .+ ar1_noise_temperature .+ temperature_0
                pulse_co2[i,:]         = sneasybrick_pulse[:ccm, :atmco2] .+ car1_noise_co2
                pulse_gmsl[i,:]        = sneasybrick_pulse[:global_sea_level, :sea_level_rise] .- mean(sneasybrick_pulse[:global_sea_level, :sea_level_rise][sealevel_norm_indices_1961_1990]) .+ ar1_noise_gmsl

            catch

                # Set values to -9999999.99 if non-physical parameter samples produce a model error.
                base_temperature[i,:]  .= -9999999.99
                base_co2[i,:]          .= -9999999.99
                base_ocean_heat[i,:]   .= -9999999.99
                base_oceanco2[i,:]     .= -9999999.99
                base_glaciers[i,:]     .= -9999999.99
                base_greenland[i,:]    .= -9999999.99
                base_antarctic[i,:]    .= -9999999.99
                base_thermal_sl[i,:]   .= -9999999.99
                base_gmsl[i,:]         .= -9999999.99

                pulse_temperature[i,:] .= -9999999.99
                pulse_co2[i,:]         .= -9999999.99
                pulse_gmsl[i,:]        .= -9999999.99
            end
        end

        # Identify model indices that caused a model error or yield non-physical outcomes (i.e. strongly negative temperatures in 2300 under RCP 8.5).
        error_indices = findall(x-> x < -10.0, base_temperature[:,end])
        good_indices  = findall(!in(error_indices), collect(1:number_samples))

        # Calculate credible intervals for base model projections that did not produce a model error
        ci_temperature = get_confidence_interval(collect(start_year:end_year), base_temperature[good_indices,:], ci_interval_1, ci_interval_2)
        ci_co2         = get_confidence_interval(collect(start_year:end_year), base_co2[good_indices,:],         ci_interval_1, ci_interval_2)
        ci_ocean_heat  = get_confidence_interval(collect(start_year:end_year), base_ocean_heat[good_indices,:],  ci_interval_1, ci_interval_2)
        ci_oceanco2    = get_confidence_interval(collect(start_year:end_year), base_oceanco2[good_indices,:],    ci_interval_1, ci_interval_2)
        ci_glaciers    = get_confidence_interval(collect(start_year:end_year), base_glaciers[good_indices,:],    ci_interval_1, ci_interval_2)
        ci_greenland   = get_confidence_interval(collect(start_year:end_year), base_greenland[good_indices,:],   ci_interval_1, ci_interval_2)
        ci_antarctic   = get_confidence_interval(collect(start_year:end_year), base_antarctic[good_indices,:],   ci_interval_1, ci_interval_2)
        ci_thermal_sl  = get_confidence_interval(collect(start_year:end_year), base_thermal_sl[good_indices,:],  ci_interval_1, ci_interval_2)
        ci_gmsl        = get_confidence_interval(collect(start_year:end_year), base_gmsl[good_indices,:],        ci_interval_1, ci_interval_2)

        # Return projections that did not cause model errors and indices for parameters that produced successful/error runs.
        return base_temperature[good_indices,:], base_co2[good_indices,:], base_ocean_heat[good_indices,:], base_oceanco2[good_indices,:],
               base_glaciers[good_indices,:], base_greenland[good_indices,:], base_antarctic[good_indices,:], base_thermal_sl[good_indices,:], base_gmsl[good_indices,:],
               pulse_temperature[good_indices,:], pulse_co2[good_indices,:], pulse_gmsl[good_indices,:],
               ci_temperature, ci_co2, ci_ocean_heat, ci_oceanco2, ci_glaciers, ci_greenland, ci_antarctic, ci_thermal_sl, ci_gmsl,
               error_indices, good_indices
    end

    # Return function with user model specifications.
    return sneasybrick_baseline
end
