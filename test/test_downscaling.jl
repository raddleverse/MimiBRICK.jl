module TestModifications

# Test use cases for the downscaling function in MimiBRICK.jl.  The `src/downscale.jl` 
# file contains routines to downscale the BRICK global sea level projections to local. 
# This uses the sea-level "fingerprints" of [Slangen et al. (2014)](https://link.springer.com/article/10.1007/s10584-014-1080-9). 

using MimiBRICK
using Test

# Lat and Lon for New York City
lat = 40.7128 # deg N
lon = 360-74.0060 # 74.0060 deg W
results_dir = joinpath(@__DIR__, "test_data/results/my_brick_results_1K_18-05-2022")

# specification 1
# testing hindcast ensemble
years, lsl_hind_ens = MimiBRICK.downscale_brick(lon=lon, lat=lat, results_dir=results_dir, proj_or_hind="hind", ensemble_or_map="ensemble", model_config="brick", rcp_scenario="RCP26")
@test size(lsl_hind_ens)[1]==168
@test size(lsl_hind_ens)[2]==100
@test length(years)==168
@test all([isa(lsl_hind_ens[i,end],Number) for i=1:length(years)])
@test all([isa(years[i],Number) for i=1:length(years)])

# testing projections ensemble
years, lsl_proj_ens = MimiBRICK.downscale_brick(lon=lon, lat=lat, results_dir=results_dir, proj_or_hind="proj", ensemble_or_map="ensemble", model_config="brick", rcp_scenario="RCP85")
@test size(lsl_proj_ens)[1]==451
@test size(lsl_proj_ens)[2]==100
@test length(years)==451
@test all([isa(lsl_proj_ens[i,end],Number) for i=1:length(years)])
@test all([isa(years[i],Number) for i=1:length(years)])

# testing projections with MAP
years, lsl_map = MimiBRICK.downscale_brick(lon=lon, lat=lat, results_dir=results_dir, proj_or_hind="proj", ensemble_or_map="map", model_config="brick", rcp_scenario="RCP85")
@test ndims(lsl_map)==1
@test length(lsl_map)==451
@test length(years)==451
@test all([isa(lsl_map[i],Number) for i=1:length(years)])
@test all([isa(years[i],Number) for i=1:length(years)])

end
