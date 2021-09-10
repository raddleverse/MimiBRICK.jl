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

## DOECLIM+BRICK
m = MimiBRICK_DOECLIM.create_brick_doeclim()
run(m)

## SNEASY+BRICK
m = create_sneasy_brick()
run(m)

##==============================================================================
