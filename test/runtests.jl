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


##==============================================================================
## With parameters read from a CSV file

## BRICK standalone
##  --> with SNEASY forcings for Temperature and Ocean heat uptake

#TODO - check what the parameter values used are - can we get these on a CSV file for reading/writing?
m = MimiBRICK.get_model()
run(m)

## DOECLIM+BRICK

#TODO - check what the parameter values used are - can we get these on a CSV file for reading/writing?
#function create_brick_doeclim(;rcp_scenario::String = "RCP85", start_year::Int=1850, end_year::Int=2300)
m = MimiBRICK_DOECLIM.create_brick_doeclim()
run(m)

## SNEASY+BRICK

#TODO - check what the parameter values used are - can we get these on a CSV file for reading/writing?
# function create_sneasy_brick(; rcp_scenario::String="RCP85", start_year::Int=1850, end_year::Int=2020)
m = create_sneasy_brick()
run(m)

##==============================================================================
