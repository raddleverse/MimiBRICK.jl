module TestDownscaling

using MimiBRICK
using Test
import MimiBRICK: get_model, create_sneasy_brick, create_brick_doeclim, downscale_brick

##==============================================================================
## Checking downscaling routine

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

end