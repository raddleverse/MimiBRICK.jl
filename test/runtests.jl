using TestItemRunner

include("test_default.jl")
include("test_modifications.jl")
include("test_calibration.jl")
include("test_downscaling.jl")

@run_package_tests
