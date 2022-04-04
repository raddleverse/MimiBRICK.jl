# Other packages
using CSV
using CSVFiles
using DataFrames
using Dates
using Distributions
using KernelDensity
using LinearAlgebra
using Mimi
using MimiBRICK
using MimiSNEASY
using NetCDF
using RobustAdaptiveMetropolisSampler
using MCMCDiagnostics
using Random
using StatsBase
using Test

srcdir = joinpath(@__DIR__, "..", "src")
include(joinpath(srcdir,"MimiBRICK_DOECLIM.jl"))
include(joinpath(srcdir,"create_models","SNEASY_BRICK.jl"))

# The land water storage is probabilistic, so global sea level change, AIS, and
# land water storage won't be exact every time. Check total GMSL against the sum
# of the components, and length of output match, as those won't be dependent on
# any future changes to model components.
testtol = 1E-6

##==============================================================================
## With default parameters

## BRICK standalone
m = MimiBRICK.get_model()
run(m)
## Run tests for the out-of-box model with default parameters and forcing
@test length(m[:global_sea_level,:sea_level_rise]) == 171
tot = m[:antarctic_icesheet,:ais_sea_level][end] + m[:glaciers_small_icecaps,:gsic_sea_level][end] +
      m[:greenland_icesheet,:greenland_sea_level][end] + m[:thermal_expansion,:te_sea_level][end] +
      m[:landwater_storage,:lws_sea_level][end]
@test m[:global_sea_level,:sea_level_rise][end] ≈ tot atol = 0.000001

## DOECLIM+BRICK
m = MimiBRICK_DOECLIM.create_brick_doeclim()
run(m)
@test length(m[:global_sea_level,:sea_level_rise]) == 171
tot = m[:antarctic_icesheet,:ais_sea_level][end] + m[:glaciers_small_icecaps,:gsic_sea_level][end] +
      m[:greenland_icesheet,:greenland_sea_level][end] + m[:thermal_expansion,:te_sea_level][end] +
      m[:landwater_storage,:lws_sea_level][end]
@test m[:global_sea_level,:sea_level_rise][end] ≈ tot atol = testtol

## SNEASY+BRICK
m = create_sneasy_brick()
run(m)
println(m[:ccm,:atmco2][end])
@test length(m[:global_sea_level,:sea_level_rise]) == 171
tot = m[:antarctic_icesheet,:ais_sea_level][end] + m[:glaciers_small_icecaps,:gsic_sea_level][end] +
      m[:greenland_icesheet,:greenland_sea_level][end] + m[:thermal_expansion,:te_sea_level][end] +
      m[:landwater_storage,:lws_sea_level][end]
@test m[:global_sea_level,:sea_level_rise][end] ≈ tot atol = testtol

##==============================================================================
## Checking other RCP scenarios, time periods, model configurations

m = MimiBRICK.get_model(rcp_scenario="RCP26", end_year=2300)
run(m)
@test length(m[:global_sea_level,:sea_level_rise]) == 451
tot = m[:antarctic_icesheet,:ais_sea_level][end] + m[:glaciers_small_icecaps,:gsic_sea_level][end] +
      m[:greenland_icesheet,:greenland_sea_level][end] + m[:thermal_expansion,:te_sea_level][end] +
      m[:landwater_storage,:lws_sea_level][end]
@test m[:global_sea_level,:sea_level_rise][end] ≈ tot atol = testtol

m = MimiBRICK.get_model(rcp_scenario="RCP45", start_year=1900, end_year=2300)
run(m)
@test length(m[:global_sea_level,:sea_level_rise]) == 401
tot = m[:antarctic_icesheet,:ais_sea_level][end] + m[:glaciers_small_icecaps,:gsic_sea_level][end] +
      m[:greenland_icesheet,:greenland_sea_level][end] + m[:thermal_expansion,:te_sea_level][end] +
      m[:landwater_storage,:lws_sea_level][end]
@test m[:global_sea_level,:sea_level_rise][end] ≈ tot atol = testtol

m = MimiBRICK.get_model(rcp_scenario="RCP60", start_year=1950, end_year=2300)
run(m)
@test length(m[:global_sea_level,:sea_level_rise]) == 351
tot = m[:antarctic_icesheet,:ais_sea_level][end] + m[:glaciers_small_icecaps,:gsic_sea_level][end] +
      m[:greenland_icesheet,:greenland_sea_level][end] + m[:thermal_expansion,:te_sea_level][end] +
      m[:landwater_storage,:lws_sea_level][end]
@test m[:global_sea_level,:sea_level_rise][end] ≈ tot atol = testtol

m = MimiBRICK.get_model(rcp_scenario="RCP85", start_year=2000, end_year=2300)
run(m)
@test length(m[:global_sea_level,:sea_level_rise]) == 301
tot = m[:antarctic_icesheet,:ais_sea_level][end] + m[:glaciers_small_icecaps,:gsic_sea_level][end] +
      m[:greenland_icesheet,:greenland_sea_level][end] + m[:thermal_expansion,:te_sea_level][end] +
      m[:landwater_storage,:lws_sea_level][end]
@test m[:global_sea_level,:sea_level_rise][end] ≈ tot atol = testtol

##==============================================================================
## Checking short calibration, that it does things and isn't just perpetually stuck

calibration_start_year = 1850
calibration_end_year   = 2017
total_chain_length     = 1000
size_subsample         = 100
threshold_gr           = 1.1

# Load calibration helper functions file.
include(joinpath("..", "calibration", "calibration_helper_functions.jl"))

# Create the log-posterior functions
include(joinpath("..", "calibration", "run_historic_models", "run_brick_historic_climate.jl"))
include(joinpath("..", "calibration", "create_log_posterior_brick.jl"))
log_posterior_brick = construct_brick_log_posterior(construct_run_brick(calibration_start_year, calibration_end_year), model_start_year=calibration_start_year, calibration_end_year=calibration_end_year, joint_antarctic_prior=false)
include(joinpath("..", "calibration", "run_historic_models", "run_doeclimbrick_historic_climate.jl"))
include(joinpath("..", "calibration", "create_log_posterior_doeclimbrick.jl"))
log_posterior_doeclimbrick = construct_doeclimbrick_log_posterior(construct_run_doeclimbrick(calibration_start_year, calibration_end_year), model_start_year=calibration_start_year, calibration_end_year=calibration_end_year, joint_antarctic_prior=false)
include(joinpath("..", "calibration", "run_historic_models", "run_sneasybrick_historic_climate.jl"))
include(joinpath("..", "calibration", "create_log_posterior_sneasybrick.jl"))
log_posterior_sneasybrick = construct_sneasybrick_log_posterior(construct_run_sneasybrick(calibration_start_year, calibration_end_year), model_start_year=calibration_start_year, calibration_end_year=calibration_end_year, joint_antarctic_prior=false)

# Do the actual calibrations
include(joinpath(@__DIR__, "..", "calibration", "calibration.jl"))

# BRICK calibration
nparameters_brick = 35
x1 = run_calibration(log_posterior_brick; model_config="brick", calibration_start_year=1850, calibration_end_year=2017,
                     total_chain_length=total_chain_length, burnin_length=0, threshold_gr=threshold_gr, num_walkers=2,
                     size_subsample=size_subsample, start_from_priors=false)
@test size(x1[1])[1]==1000
@test size(x1[1])[2]==nparameters_brick
@test size(x1[5])[1]==100
@test size(x1[5])[2]==nparameters_brick
@test all([isa(x1[1][end,i],Number) for i=1:size(x1[1])[2]])
@test !all([diff(x1[1][:,1])[i] == 0 for i=1:size(x1[1])[1]-1])

# DOECLIM-BRICK calibration
nparameters_doeclimbrick = 44
x2 = run_calibration(log_posterior_doeclimbrick; model_config="doeclimbrick", calibration_start_year=1850, calibration_end_year=2017,
                     total_chain_length=total_chain_length, burnin_length=0, threshold_gr=threshold_gr, num_walkers=2,
                     size_subsample=size_subsample, start_from_priors=false)
@test size(x2[1])[1]==1000
@test size(x2[1])[2]==nparameters_doeclimbrick
@test size(x2[5])[1]==100
@test size(x2[5])[2]==nparameters_doeclimbrick
@test all([isa(x2[1][end,i],Number) for i=1:size(x2[1])[2]])
@test !all([diff(x2[1][:,1])[i] == 0 for i=1:size(x2[1])[1]-1])

# SNEASY-BRICK calibration
nparameters_sneasybrick = 51
x3 = run_calibration(log_posterior_sneasybrick; model_config="sneasybrick", calibration_start_year=1850, calibration_end_year=2017,
                     total_chain_length=total_chain_length, burnin_length=0, threshold_gr=threshold_gr, num_walkers=2,
                     size_subsample=size_subsample, start_from_priors=false)
@test size(x3[1])[1]==1000
@test size(x3[1])[2]==nparameters_sneasybrick
@test size(x3[5])[1]==100
@test size(x3[5])[2]==nparameters_sneasybrick
@test all([isa(x3[1][end,i],Number) for i=1:size(x3[1])[2]])
@test !all([diff(x3[1][:,1])[i] == 0 for i=1:size(x3[1])[1]-1])

##==============================================================================
## Checking downscaling routine

include(joinpath(@__DIR__, "..", "localslr", "downscale.jl"))

# for testing - New York City
lat = 40.7128 # deg N
lon = 360-74.0060 # 74.0060 deg W

# testing hindcast ensemble
years, lsl_hind_ens = downscale_brick(lon=lon, lat=lat, proj_or_hind="hind", ensemble_or_map="ensemble", model_config="brick", rcp_scenario="RCP26")
@test size(lsl_hind_ens)[1]==168
@test size(lsl_hind_ens)[2]==10000
@test length(years)==168
@test all([isa(lsl_hind_ens[end,i],Number) for i=1:length(years)])
@test all([isa(years[i],Number) for i=1:length(years)])

# testing projections ensemble
years, lsl_proj_ens = downscale_brick(lon=lon, lat=lat, proj_or_hind="proj", ensemble_or_map="ensemble", model_config="brick", rcp_scenario="RCP85")
@test size(lsl_proj_ens)[1]==451
@test size(lsl_proj_ens)[2]==10000
@test length(years)==451
@test all([isa(lsl_proj_ens[end,i],Number) for i=1:length(years)])
@test all([isa(years[i],Number) for i=1:length(years)])

# testing projections with MAP
years, lsl_map = downscale_brick(lon=lon, lat=lat, proj_or_hind="proj", ensemble_or_map="map", model_config="brick", rcp_scenario="RCP60")
@test ndims(lsl_map)==1
@test length(lsl_map)==451
@test length(years)==451
@test all([isa(lsl_map[i],Number) for i=1:length(years)])
@test all([isa(years[i],Number) for i=1:length(years)])

##==============================================================================
