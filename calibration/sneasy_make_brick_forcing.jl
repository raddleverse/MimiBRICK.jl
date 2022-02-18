# Activate the project for the paper and make sure all packages we need are installed.
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

# Other packages
using CSV
using DataFrames
using MimiBRICK
using MimiSNEASY

srcdir = joinpath(@__DIR__, "..", "src")
include(joinpath(srcdir,"MimiBRICK_DOECLIM.jl"))
include(joinpath(srcdir,"create_models","SNEASY_BRICK.jl"))

##==============================================================================
## Read output from SNEASY RCP projections

#TODO

## Compute mean temperature and ocean heat responses

#TODO

## Normalize temperature relative to 1850-1870

#TODO

## Write CSV files

#TODO


##==============================================================================
