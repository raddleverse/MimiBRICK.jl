##==============================================================================
## Script for running BRICK (standalone) over the historic period and save the
## hindcast results to CSV files.
##==============================================================================


##==============================================================================
## Initial set-up

# Activate the project for the paper and make sure all packages we need are installed.
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

# Other packages
using Test
using CSV
using DataFrames
using MimiBRICK
using MimiSNEASY
using LinearAlgebra

srcdir = joinpath(@__DIR__, "..", "src")
include(joinpath(srcdir,"MimiBRICK_DOECLIM.jl"))
include(joinpath(srcdir,"create_models","SNEASY_BRICK.jl"))
include("calibration_helper_functions.jl")
include(joinpath(@__DIR__, "..", "calibration", "helper_functions.jl"))

# Model configuration
# --> Possible options: (1) "brick", (2) "doeclimbrick", (3) "sneasybrick"
model_config = "brick"

# Set years for model calibration.
start_year = calibration_start_year = 1850
end_year = calibration_end_year = 2017
model_years  = collect(calibration_start_year:calibration_end_year)
number_years = length(model_years)

##==============================================================================
## Set paths for results files - subsample of model parameters

filename_brick_parameters = joinpath(@__DIR__, "..", "results", "my_brick_results_20M_09-02-2022", "parameters_subsample.csv")
filename_doeclimbrick_parameters = joinpath(@__DIR__, "..", "results", "my_doeclimbrick_results_20M_09-02-2022", "parameters_subsample.csv")doeclim
filename_sneasybrick_parameters = joinpath(@__DIR__, "..", "results", "my_sneasybrick_results_20M_09-02-2022", "parameters_subsample.csv")

##==============================================================================
## Read subsample of parameters

if model_config == "brick"
    parameters = DataFrame(load(filename_brick_parameters))
elseif model_config == "doeclimbrick"
    parameters = DataFrame(load(filename_doeclimbrick_parameters))
elseif model_config == "sneasybrick"
    parameters = DataFrame(load(filename_sneasybrick_parameters))
end
num_ens = size(parameters)[1]
num_par = size(parameters)[2]
parnames = names(parameters)

##==============================================================================
## Run model at each set

# Get model instance
m = MimiBRICK.get_model()

# Initialize arrays to save the model components

# Pre-allocate arrays to store (SNEASY/DOECLIM+)BRICK results.
glaciers     = zeros(Union{Missing, Float64}, num_ens, number_years)
greenland    = zeros(Union{Missing, Float64}, num_ens, number_years)
thermal_sl   = zeros(Union{Missing, Float64}, num_ens, number_years)
antarctic    = zeros(Union{Missing, Float64}, num_ens, number_years)
gmsl         = zeros(Union{Missing, Float64}, num_ens, number_years)
# Pre-allocate vectors to hold simulated CAR(1) & AR(1) with measurement error noise.
ar1_noise_glaciers    = zeros(number_years)
ar1_noise_greenland   = zeros(number_years)
ar1_noise_antarctic   = zeros(number_years)
ar1_noise_gmsl        = zeros(number_years)
# Replicate errors for years without observations over model time horizon (used for simulating AR1 noise).
obs_error_glaciers    = replicate_errors(start_year, end_year, calibration_data.glaciers_sigma)
obs_error_greenland   = replicate_errors(start_year, end_year, calibration_data.merged_greenland_sigma)
obs_error_antarctic   = replicate_errors(start_year, end_year, calibration_data.antarctic_imbie_sigma)
obs_error_gmsl        = replicate_errors(start_year, end_year, calibration_data.gmsl_sigma)
if (model_config == "doeclimbrick") | (model_config == "sneasybrick")
    temperature  = zeros(Union{Missing, Float64}, num_ens, number_years)
    ocean_heat   = zeros(Union{Missing, Float64}, num_ens, number_years)
    ar1_noise_temperature = zeros(number_years)
    ar1_noise_ocean_heat  = zeros(number_years)
    obs_error_temperature = replicate_errors(start_year, end_year, calibration_data.hadcrut_temperature_sigma)
    obs_error_oceanheat   = replicate_errors(start_year, end_year, calibration_data.ocean_heat_sigma)
    if (model_config == "sneasybrick")
        co2          = zeros(Union{Missing, Float64}, num_ens, number_years)
        oceanco2     = zeros(Union{Missing, Float64}, num_ens, number_years)
        normal_noise_oceanco2 = zeros(number_years)
        car1_noise_co2        = zeros(number_years)
    end
end

# Get indices needed to normalize temperature anomalies relative to 1861-1880 mean.
temperature_norm_indices = findall((in)(1861:1880), calibration_start_year:calibration_end_year)

# Get indices needed to normalize all sea level rise sources.
sealevel_norm_indices_1961_1990 = findall((in)(1961:1990), calibration_start_year:calibration_end_year)
sealevel_norm_indices_1992_2001 = findall((in)(1992:2001), calibration_start_year:calibration_end_year)

# Load calibration data from 1765-2017 (measurement errors used in simulated noise).
calibration_data, obs_antarctic_trends, obs_thermal_trends = load_calibration_data(start_year, 2017)

# Loop over parameters and run it

#TODO
for i = 1:num_ens

    # TODO - set all parameters
    # Start with the BRICK ones, common to all 3 configurations (BRICK, DOECLIM-BRICK, SNEASY-BRICK)
    # ----- Antarctic Ocean ----- #
    update_param!(m, :antarctic_ocean, :anto_α, parameters[i, findall(x->x=="anto_alpha",parnames)][1])
    update_param!(m, :antarctic_ocean, :anto_β, parameters[i, findall(x->x=="anto_beta",parnames)][1])

    # ----- Antarctic Ice Sheet ----- #
    update_param!(m, :antarctic_icesheet, :ais_sea_level₀, parameters[i, findall(x->x=="antarctic_s0",parnames)][1])
    update_param!(m, :antarctic_icesheet, :ais_bedheight₀, parameters[i, findall(x->x=="antarctic_bed_height0",parnames)][1])
    update_param!(m, :antarctic_icesheet, :ais_slope, parameters[i, findall(x->x=="antarctic_slope",parnames)][1])
    update_param!(m, :antarctic_icesheet, :ais_μ, parameters[i, findall(x->x=="antarctic_mu",parnames)][1])
    update_param!(m, :antarctic_icesheet, :ais_runoffline_snowheight₀, parameters[i, findall(x->x=="antarctic_runoff_height0",parnames)][1])
    update_param!(m, :antarctic_icesheet, :ais_c, parameters[i, findall(x->x=="antarctic_c",parnames)][1])
    update_param!(m, :antarctic_icesheet, :ais_precipitation₀, parameters[i, findall(x->x=="antarctic_precip0",parnames)][1])
    update_param!(m, :antarctic_icesheet, :ais_κ, parameters[i, findall(x->x=="antarctic_kappa",parnames)][1])
    update_param!(m, :antarctic_icesheet, :ais_ν, parameters[i, findall(x->x=="antarctic_nu",parnames)][1])
    update_param!(m, :antarctic_icesheet, :ais_iceflow₀, parameters[i, findall(x->x=="antarctic_flow0",parnames)][1])
    update_param!(m, :antarctic_icesheet, :ais_γ, parameters[i, findall(x->x=="antarctic_gamma",parnames)][1])
    update_param!(m, :antarctic_icesheet, :ais_α, parameters[i, findall(x->x=="antarctic_alpha",parnames)][1])
    update_param!(m, :antarctic_icesheet, :temperature_threshold, parameters[i, findall(x->x=="antarctic_temp_threshold",parnames)][1])
    update_param!(m, :antarctic_icesheet, :λ, parameters[i, findall(x->x=="antarctic_lambda",parnames)][1])

    # ----- Glaciers & Small Ice Caps ----- #
    update_param!(m, :glaciers_small_icecaps, :gsic_β₀, parameters[i, findall(x->x=="glaciers_beta0",parnames)][1])
    update_param!(m, :glaciers_small_icecaps, :gsic_v₀, parameters[i, findall(x->x=="glaciers_v0",parnames)][1])
    update_param!(m, :glaciers_small_icecaps, :gsic_s₀, parameters[i, findall(x->x=="glaciers_s0",parnames)][1])
    update_param!(m, :glaciers_small_icecaps, :gsic_n, parameters[i, findall(x->x=="glaciers_n",parnames)][1])
    #update_param!(m, :glaciers_small_icecaps, :gsic_teq, parameters[i, findall(x->x=="anto_beta",parnames)][1])

    # ----- Greenland Ice Sheet ----- #
    update_param!(m, :greenland_icesheet, :greenland_a, parameters[i, findall(x->x=="greenland_a",parnames)][1])
    update_param!(m, :greenland_icesheet, :greenland_b, parameters[i, findall(x->x=="greenland_b",parnames)][1])
    update_param!(m, :greenland_icesheet, :greenland_α, parameters[i, findall(x->x=="greenland_alpha",parnames)][1])
    update_param!(m, :greenland_icesheet, :greenland_β, parameters[i, findall(x->x=="greenland_beta",parnames)][1])
    update_param!(m, :greenland_icesheet, :greenland_v₀, parameters[i, findall(x->x=="greenland_v0",parnames)][1])

    # ----- Thermal Expansion ----- #
    update_param!(m, :thermal_expansion, :te_α, parameters[i, findall(x->x=="thermal_alpha",parnames)][1])
    update_param!(m, :thermal_expansion, :te_s₀, parameters[i, findall(x->x=="thermal_s0",parnames)][1])

    if (model_config == "doeclimbrick") | (model_config == "sneasybrick")
        # add the DOECLIM/SNEASY common parameters
        # TODO
        if (model_config == "sneasybrick")
            # add the SNEASY-only parameters
            # TODO
        end
    end

    # run model
    run(m)

    # AR1 noise
    # TODO - need to modify these
    if false
        # Create noise to superimpose on results using calibrated statistical parameters and measurement noise (note: both models use same estimated noise for each sample).
        σ_glaciers               = parameters[i, findall(x->x=="sd_glaciers",parnames)][1]
        σ_greenland              = parameters[i, findall(x->x=="sd_greenland",parnames)][1]
        σ_antarctic              = parameters[i, findall(x->x=="sd_antarctic",parnames)][1]
        σ_gmsl                   = parameters[i, findall(x->x=="sd_gmsl",parnames)][1]
        ρ_glaciers               = parameters[i, findall(x->x=="rho_glaciers",parnames)][1]
        ρ_greenland              = parameters[i, findall(x->x=="rho_glaciers",parnames)][1]
        ρ_antarctic              = parameters[i, findall(x->x=="rho_glaciers",parnames)][1]
        ρ_gmsl                   = parameters[i, findall(x->x=="rho_glaciers",parnames)][1]
        ar1_noise_glaciers[:]    = simulate_ar1_noise(number_years, σ_glaciers,    ρ_glaciers,    obs_error_glaciers)
        ar1_noise_greenland[:]   = simulate_ar1_noise(number_years, σ_greenland,   ρ_greenland,   obs_error_greenland)
        ar1_noise_antarctic[:]   = simulate_ar1_noise(number_years, σ_antarctic,   ρ_antarctic,   obs_error_antarctic)
        ar1_noise_gmsl[:]        = simulate_ar1_noise(number_years, σ_gmsl,        ρ_gmsl,        obs_error_gmsl)
        if (model_config == "doeclimbrick") | (model_config == "sneasybrick")
            σ_temperature            = parameters[i, findall(x->x=="sd_temp",parnames)][1]
            σ_ocean_heat             = parameters[i, findall(x->x=="sd_ocean_heat",parnames)][1]
            ρ_temperature            = parameters[i, findall(x->x=="rho_temperature",parnames)][1]
            ρ_ocean_heat             = parameters[i, findall(x->x=="rho_ocean_heat",parnames)][1]
            ar1_noise_temperature[:] = simulate_ar1_noise(number_years, σ_temperature, ρ_temperature, obs_error_temperature)
            ar1_noise_ocean_heat[:]  = simulate_ar1_noise(number_years, σ_ocean_heat,  ρ_ocean_heat,  obs_error_oceanheat)
            if (model_config == "sneasybrick")
                σ²_white_noise_CO₂       = parameters[i, findall(x->x=="sigma_whitenoise_co2",parnames)][1]
                normal_noise_oceanco2[:] = rand(Normal(0,0.4*sqrt(10)), number_years)
                # CO₂ uses CAR(1) statistical process parameters.
                car1_noise_co2[:] = simulate_car1_noise(number_years, α₀_CO₂, σ²_white_noise_CO₂, obs_error_co2)
            end
        end

        # Normalize relative to appropriate time period, and superimpose statistical noise where appropriate

# HERE NOW - SORTING OUT THE SIZE MISMATCH
    # HERE NOW - SORTING OUT THE SIZE MISMATCH
        # HERE NOW - SORTING OUT THE SIZE MISMATCH
# HERE NOW - SORTING OUT THE SIZE MISMATCH
    # HERE NOW - SORTING OUT THE SIZE MISMATCH
        # HERE NOW - SORTING OUT THE SIZE MISMATCH

        glaciers[i,:]    = m[:glaciers_small_icecaps, :gsic_sea_level] .- mean(m[:glaciers_small_icecaps, :gsic_sea_level][sealevel_norm_indices_1961_1990]) .+ ar1_noise_glaciers
        greenland[i,:]   = m[:greenland_icesheet, :greenland_sea_level] .- mean(m[:greenland_icesheet, :greenland_sea_level][sealevel_norm_indices_1992_2001]) .+ ar1_noise_greenland
        antarctic[i,:]   = m[:antarctic_icesheet, :ais_sea_level] .- mean(m[:antarctic_icesheet, :ais_sea_level][sealevel_norm_indices_1992_2001]) .+ ar1_noise_antarctic
        thermal_sl[i,:]  = m[:thermal_expansion, :te_sea_level]
        gmsl[i,:]        = m[:global_sea_level, :sea_level_rise] .- mean(m[:global_sea_level, :sea_level_rise][sealevel_norm_indices_1961_1990]) .+ ar1_noise_gmsl
        if (model_config == "doeclimbrick") | (model_config == "sneasybrick")
            temperature[i,:] = m[:doeclim, :temp] .- mean(m[:doeclim, :temp][temperature_norm_indices]) .+ ar1_noise_temperature .+ temperature_0
            ocean_heat[i,:]  = m[:doeclim, :heat_mixed] .+ m[:doeclim, :heat_interior] .+ ar1_noise_ocean_heat .+ ocean_heat_0
            if (model_config == "sneasybrick")
                co2[i,:]         = m[:ccm, :atmco2] .+ car1_noise_co2
                oceanco2[i,:]    = m[:ccm, :atm_oc_flux] .+ normal_noise_oceanco2
            end
        end
    end


    # save results
end


##==============================================================================
## Save output

# create results file names based on the initial parameter data set

#TODO (edit these - just placeholders from the CIAM code)
outtrialsname = joinpath(postprocessing_outputdir, "trials_$(runname).csv")
CSV.write(outtrialsname, outtrials)




##==============================================================================
## End
##==============================================================================
