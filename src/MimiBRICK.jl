module MimiBRICK

# Load required packages
using CSVFiles, DataFrames, Distributions, Mimi

# Load BRICK Mimi components.
include(joinpath("components", "antarctic_icesheet_component.jl"))
include(joinpath("components", "antarctic_ocean_component.jl"))
include(joinpath("components", "glaciers_small_icecaps_component.jl"))
include(joinpath("components", "global_sea_level_component.jl"))
include(joinpath("components", "greenland_icesheet_component.jl"))
include(joinpath("components", "landwater_storage_component.jl"))
include(joinpath("components", "thermal_expansion_component.jl"))

# Load creation functions for two other model variants
include(joinpath("create_models", "SNEASY_BRICK.jl"))
include(joinpath("create_models", "BRICK_DOECLIM.jl"))

# Load other helper functions
include(joinpath("downscale.jl"))

# Load calibration functions
include(joinpath("calibration/main_calibration.jl"))

# Load get_model function
include("get_model.jl")

end # Module
