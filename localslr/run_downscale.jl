
# Activate the project for the paper and make sure all packages we need are installed.
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include(joinpath(@__DIR__, "..", "localslr", "downscale.jl"))

# Example - New York City
lat = 40.7128 # deg N
lon = 360-74.0060 # 74.0060 deg W
model_config = "brick"
proj_or_hind = "proj"
rcp_scenario = "RCP85"
ensemble_or_map = "map"

# Example: using RCP8.5 projections and the maximum a posteriori ensemble member under
years, lsl = downscale_brick(lon=lon, lat=lat, proj_or_hind=proj_or_hind, ensemble_or_map=ensemble_or_map, model_config=model_config, rcp_scenario=rcp_scenario)
