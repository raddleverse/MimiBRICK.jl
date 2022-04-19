using Test
using MimiBRICK

# @testset "MimiBRICK" begin 

      @info("test_default.jl")
      @time include("test_default.jl")

      @info("test_modifications.jl")
      @time include("test_modifications.jl")

      @info("test_calibration.jl")
      @time include("test_calibration.jl")

      @info("test_downscaling.jl")
      # @time include("test_downscaling.jl") # TODO generalize this testing so it doesn't depend on local files
# end
