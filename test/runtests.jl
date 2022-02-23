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

srcdir = joinpath(@__DIR__, "..", "src")
include(joinpath(srcdir,"MimiBRICK_DOECLIM.jl"))
include(joinpath(srcdir,"create_models","SNEASY_BRICK.jl"))

##==============================================================================
## With default parameters

## BRICK standalone
m = MimiBRICK.get_model()
run(m)
## Run tests for the out-of-box model with default parameters and forcing
@test length(m[:global_sea_level,:sea_level_rise]) == 451
@test m[:global_sea_level,:sea_level_rise][end] == 8.83141199959874
@test m[:antarctic_icesheet,:ais_sea_level][end] == 2.565599697283293
@test m[:glaciers_small_icecaps,:gsic_sea_level][end] == 0.37599982678123345
@test m[:greenland_icesheet,:greenland_sea_level][end] == 4.617681707903225
@test m[:thermal_expansion,:te_sea_level][end] == 1.1862520569095971

## DOECLIM+BRICK
m = MimiBRICK_DOECLIM.create_brick_doeclim()
run(m)

## SNEASY+BRICK
m = create_sneasy_brick()
run(m)

##==============================================================================
