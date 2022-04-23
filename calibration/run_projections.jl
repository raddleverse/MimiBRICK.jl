##==============================================================================
## Script for running BRICK (standalone, or with DOECLIM or SNEASY) over the
## projections period and save the results to CSV files.
## Projections under RCP scenarios (2.6, 4.5, 6.0, or 8.5)
##==============================================================================


##==============================================================================
## Initial set-up

# Other packages
using Test
using CSV
using DataFrames
using MimiSNEASY
using LinearAlgebra

include(joinpath("calibration_helper_functions.jl"))
include(joinpath("..", "calibration", "helper_functions.jl"))
outdir = joinpath(@__DIR__, "..", "results")

# Model configuration
# --> Possible options: (1) "brick", (2) "doeclimbrick", (3) "sneasybrick"
model_config = "sneasybrick"

# RCP scenario
# --> Possible options: (1) RCP26, (2) RCP45, (3) RCP60, (4) RCP85
rcp_scenario = "RCP85"

##==============================================================================
## Should not need to mess around with anything below here
##==============================================================================

# Set years for model calibration.
start_year = 1850 # BRICK limited to 1850-2300
end_year = 2300
model_years  = collect(start_year:end_year)
num_years = length(model_years)


##==============================================================================
## Set paths for results files - subsample of model parameters, and associated log-posterior scores

# for BRICK
dir_brick = joinpath(@__DIR__, "..", "results", "my_brick_results_20M_20-02-2022")
# for DOECLIM-BRICK
dir_doeclimbrick = joinpath(@__DIR__, "..", "results", "my_doeclimbrick_results_20M_19-02-2022")
# for SNEASY-BRICK
dir_sneasybrick = joinpath(@__DIR__, "..", "results", "my_sneasybrick_results_20M_19-02-2022")

##==============================================================================
## Modify below here at your own risk
##==============================================================================


##==============================================================================
## Read subsample of parameters

if model_config == "brick"
    dir_output = dir_brick
elseif model_config == "doeclimbrick"
    dir_output = dir_doeclimbrick
elseif model_config == "sneasybrick"
    dir_output = dir_sneasybrick
end
appen = "$(model_config)_$(dir_output[(findfirst("results_", dir_output)[end]+1):length(dir_output)])"
filename_parameters = joinpath(dir_output, "parameters_subsample_$(appen).csv")
filename_logpost    = joinpath(dir_output, "log_post_subsample_$(appen).csv")
parameters = DataFrame(load(filename_parameters))
logpost = DataFrame(load(filename_logpost))[!,:log_post]
num_ens = size(parameters)[1]
num_par = size(parameters)[2]
parnames = names(parameters)

##==============================================================================
## Run model at each set

# Get model instance
if model_config=="brick"
    m = MimiBRICK.get_model(rcp_scenario=rcp_scenario, start_year=start_year, end_year=end_year)
elseif model_config=="doeclimbrick"
    m = MimiBRICK.create_brick_doeclim(rcp_scenario=rcp_scenario, start_year=start_year, end_year=end_year)
elseif model_config=="sneasybrick"
    m = MimiBRICK.create_sneasy_brick(rcp_scenario=rcp_scenario, start_year=start_year, end_year=end_year)
end

# Load calibration data from 1765-2017 (measurement errors used in simulated noise).
calibration_data, obs_antarctic_trends, obs_thermal_trends = load_calibration_data(start_year, 2017)

# Initialize arrays to save the model components

# Pre-allocate arrays to store (SNEASY/DOECLIM+)BRICK results.
glaciers     = zeros(Union{Missing, Float64}, num_ens, num_years)
greenland    = zeros(Union{Missing, Float64}, num_ens, num_years)
thermal_sl   = zeros(Union{Missing, Float64}, num_ens, num_years)
antarctic    = zeros(Union{Missing, Float64}, num_ens, num_years)
gmsl         = zeros(Union{Missing, Float64}, num_ens, num_years)
# Also need to calculate landwater storage contribution so it is the same between base and pulse runs.
landwater_storage_sl = zeros(Union{Missing, Float64}, num_ens, num_years)
# Pre-allocate vectors to hold simulated CAR(1) & AR(1) with measurement error noise.
ar1_noise_glaciers    = zeros(num_years)
ar1_noise_greenland   = zeros(num_years)
ar1_noise_antarctic   = zeros(num_years)
ar1_noise_gmsl        = zeros(num_years)
# Replicate errors for years without observations over model time horizon (used for simulating AR1 noise).
obs_error_glaciers    = replicate_errors(start_year, end_year, calibration_data.glaciers_sigma)
obs_error_greenland   = replicate_errors(start_year, end_year, calibration_data.merged_greenland_sigma)
obs_error_antarctic   = replicate_errors(start_year, end_year, calibration_data.antarctic_imbie_sigma)
obs_error_gmsl        = replicate_errors(start_year, end_year, calibration_data.gmsl_sigma)
if (model_config == "doeclimbrick") | (model_config == "sneasybrick")
    temperature  = zeros(Union{Missing, Float64}, num_ens, num_years)
    ocean_heat   = zeros(Union{Missing, Float64}, num_ens, num_years)
    ar1_noise_temperature = zeros(num_years)
    ar1_noise_ocean_heat  = zeros(num_years)
    obs_error_temperature = replicate_errors(start_year, end_year, calibration_data.hadcrut_temperature_sigma)
    obs_error_oceanheat   = replicate_errors(start_year, end_year, calibration_data.ocean_heat_sigma)
    if (model_config == "sneasybrick")
        co2          = zeros(Union{Missing, Float64}, num_ens, num_years)
        oceanco2     = zeros(Union{Missing, Float64}, num_ens, num_years)
        normal_noise_oceanco2 = zeros(num_years)
        car1_noise_co2        = zeros(num_years)
        # For CAR(1) noise, set constant CO₂ observations errors for start year-1958 (Law Dome), and 1959-end year (Mauna Loa).
        obs_error_co2 = ones(num_years) .* unique(skipmissing(calibration_data.maunaloa_co2_sigma))[1]
        obs_error_co2[1:findfirst(isequal(1958), model_years)] .= unique(skipmissing(calibration_data.lawdome_co2_sigma))[1]
    end
end

# Get indices needed to normalize temperature anomalies relative to 1861-1880 mean.
temperature_norm_indices = findall((in)(1861:1880), start_year:end_year)

# Get indices needed to normalize all sea level rise sources.
sealevel_norm_indices_1961_1990 = findall((in)(1961:1990), start_year:end_year)
sealevel_norm_indices_1992_2001 = findall((in)(1992:2001), start_year:end_year)

# Loop over parameters and run the model
for i = 1:num_ens

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

    # Calculate land water storage contribution to sea level rise (sampled from Normal distribution) and set same scenario for base and pulse runs.
    landwater_storage_sl[i,:] = rand(Normal(0.0003, 0.00018), num_years)
    update_param!(m, :landwater_storage, :lws_random_sample, landwater_storage_sl[i,:])
    update_param!(m, :landwater_storage, :lws_random_sample, landwater_storage_sl[i,:])

    if (model_config == "doeclimbrick") | (model_config == "sneasybrick")
        # add the DOECLIM/SNEASY common parameters
        update_param!(m, :doeclim, :t2co, parameters[i, findall(x->x=="climate_sensitivity",parnames)][1])
        update_param!(m, :doeclim, :kappa, parameters[i, findall(x->x=="heat_diffusivity",parnames)][1])
        update_param!(m, :radiativeforcing, :alpha, parameters[i, findall(x->x=="rf_scale_aerosol",parnames)][1])
        if (model_config == "sneasybrick")
            # add the SNEASY-only parameters
            update_param!(m, :model_CO₂_0, parameters[i, findall(x->x=="CO2_0",parnames)][1])
            update_param!(m, :ccm, :Q10, parameters[i, findall(x->x=="Q10",parnames)][1])
            update_param!(m, :ccm, :Beta, parameters[i, findall(x->x=="CO2_fertilization",parnames)][1])
            update_param!(m, :ccm, :Eta, parameters[i, findall(x->x=="CO2_diffusivity",parnames)][1])
            update_param!(m, :rfco2, :N₂O_0, parameters[i, findall(x->x=="N2O_0",parnames)][1])
        end
    end

    # run model
    run(m)

    # Statistical noise models
    # Create noise to superimpose on results using calibrated statistical parameters and measurement noise (note: both models use same estimated noise for each sample).
    σ_glaciers               = parameters[i, findall(x->x=="sd_glaciers",parnames)][1]
    σ_greenland              = parameters[i, findall(x->x=="sd_greenland",parnames)][1]
    σ_antarctic              = parameters[i, findall(x->x=="sd_antarctic",parnames)][1]
    σ_gmsl                   = parameters[i, findall(x->x=="sd_gmsl",parnames)][1]
    ρ_glaciers               = parameters[i, findall(x->x=="rho_glaciers",parnames)][1]
    ρ_greenland              = parameters[i, findall(x->x=="rho_glaciers",parnames)][1]
    ρ_antarctic              = parameters[i, findall(x->x=="rho_glaciers",parnames)][1]
    ρ_gmsl                   = parameters[i, findall(x->x=="rho_glaciers",parnames)][1]
    ar1_noise_glaciers[:]    = simulate_ar1_noise(num_years, σ_glaciers,    ρ_glaciers,    obs_error_glaciers)
    ar1_noise_greenland[:]   = simulate_ar1_noise(num_years, σ_greenland,   ρ_greenland,   obs_error_greenland)
    ar1_noise_antarctic[:]   = simulate_ar1_noise(num_years, σ_antarctic,   ρ_antarctic,   obs_error_antarctic)
    ar1_noise_gmsl[:]        = simulate_ar1_noise(num_years, σ_gmsl,        ρ_gmsl,        obs_error_gmsl)
    if (model_config == "doeclimbrick") | (model_config == "sneasybrick")
        temperature_0            = parameters[i, findall(x->x=="temperature_0",parnames)][1]
        ocean_heat_0             = parameters[i, findall(x->x=="ocean_heat_0",parnames)][1]
        σ_temperature            = parameters[i, findall(x->x=="sd_temp",parnames)][1]
        σ_ocean_heat             = parameters[i, findall(x->x=="sd_ocean_heat",parnames)][1]
        ρ_temperature            = parameters[i, findall(x->x=="rho_temperature",parnames)][1]
        ρ_ocean_heat             = parameters[i, findall(x->x=="rho_ocean_heat",parnames)][1]
        ar1_noise_temperature[:] = simulate_ar1_noise(num_years, σ_temperature, ρ_temperature, obs_error_temperature)
        ar1_noise_ocean_heat[:]  = simulate_ar1_noise(num_years, σ_ocean_heat,  ρ_ocean_heat,  obs_error_oceanheat)
        if (model_config == "sneasybrick")
            α₀_CO₂                   = parameters[i, findall(x->x=="alpha0_CO2",parnames)][1]
            σ²_white_noise_CO₂       = parameters[i, findall(x->x=="sigma_whitenoise_co2",parnames)][1]
            normal_noise_oceanco2[:] = rand(Normal(0,0.4*sqrt(10)), num_years)
            # CO₂ uses CAR(1) statistical process parameters.
            car1_noise_co2[:] = simulate_car1_noise(num_years, α₀_CO₂, σ²_white_noise_CO₂, obs_error_co2)
        end
    end
    # Normalize relative to appropriate time period, and superimpose statistical noise where appropriate
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


##==============================================================================
## Save output

# get just the specific run name that matches the directory created for the calibration output
model_tag = filename_parameters[(findfirst("results/", filename_parameters)[end]+1):(findfirst("/parameters_subsample", filename_parameters)[1]-1)]
outdir = filename_parameters[1:(findfirst("/parameters_subsample", filename_parameters)[1]-1)]
# make appropriate directory if needed
filepath_output = joinpath(outdir, "projections_csv",rcp_scenario)
mkpath(filepath_output)

# Transposing so each column is a different ensemble member, and each row is a different year
function write_output_table(field_name, appen, rcp, field_array, output_path)
    filename_output = joinpath(output_path,"projections_$(field_name)_$(rcp)_$(appen).csv")
    CSV.write(filename_output, DataFrame(field_array', :auto))
end

# Writing output tables.
write_output_table("gmsl", appen, rcp_scenario, gmsl, filepath_output)
write_output_table("landwater_storage_sl", appen, rcp_scenario, landwater_storage_sl, filepath_output)
write_output_table("glaciers", appen, rcp_scenario, glaciers, filepath_output)
write_output_table("greenland", appen, rcp_scenario, greenland, filepath_output)
write_output_table("antarctic", appen, rcp_scenario, antarctic, filepath_output)
write_output_table("thermal", appen, rcp_scenario, thermal_sl, filepath_output)
if (model_config == "doeclimbrick") | (model_config == "sneasybrick")
    write_output_table("temperature", appen, rcp_scenario, temperature, filepath_output)
    write_output_table("ocean_heat", appen, rcp_scenario, ocean_heat, filepath_output)
    if (model_config == "sneasybrick")
        write_output_table("co2", appen, rcp_scenario, co2, filepath_output)
        write_output_table("oceanco2", appen, rcp_scenario, oceanco2, filepath_output)
    end
end

# write maximum a posteriori ensemble member
idx_max = findmax(logpost)[2]
if model_config=="brick"
    colnames_out = ["YEAR","GMSL","LWS","GLAC","GIS","AIS","TE"]
elseif model_config=="doeclimbrick"
    colnames_out = ["YEAR","GMSL","LWS","GLAC","GIS","AIS","TE","TEMP","OCHEAT"]
elseif model_config=="sneasybrick"
    colnames_out = ["YEAR","GMSL","LWS","GLAC","GIS","AIS","TE","TEMP","OCHEAT","CO2","OCEANCO2"]
end
num_outputs = size(colnames_out)[1]
map_outputs = zeros(Union{Missing, Float64}, num_outputs, num_years)
map_outputs[1,:] = model_years
map_outputs[2,:] = gmsl[idx_max,:]
map_outputs[3,:] = landwater_storage_sl[idx_max,:]
map_outputs[4,:] = glaciers[idx_max,:]
map_outputs[5,:] = greenland[idx_max,:]
map_outputs[6,:] = antarctic[idx_max,:]
map_outputs[7,:] = thermal_sl[idx_max,:]
if (model_config=="doeclimbrick") | (model_config=="sneasybrick")
    map_outputs[8,:] = temperature[idx_max,:]
    map_outputs[9,:] = ocean_heat[idx_max,:]
    if (model_config=="sneasybrick")
        map_outputs[10,:] = co2[idx_max,:]
        map_outputs[11,:] = oceanco2[idx_max,:]
    end
end
df_map_outputs = DataFrame(map_outputs', colnames_out)
filename_map_outputs = joinpath(filepath_output,"projections_MAP_$(rcp_scenario)_$(appen).csv")
CSV.write(filename_map_outputs, df_map_outputs)


##==============================================================================
## End
##==============================================================================
