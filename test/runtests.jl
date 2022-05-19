using Test
using MimiBRICK

@testset "MimiBRICK" begin 

      @info("test_default.jl")
      @time include("test_default.jl")

      @info("test_modifications.jl")
      @time include("test_modifications.jl")

      @info("test_calibration.jl")
      @time include("test_calibration.jl")

end
