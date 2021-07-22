using Test
using CSV
using DataFrames
using MimiBRICK
using MimiSNEASY

include("MimiBRICK_DOECLIM.jl")
include("create_models/SNEASY_BRICK.jl")

##==============================================================================
## BRICK standalone
##  --> with SNEASY forcings for Temperature and Ocean heat uptake

#TODO - check what the parameter values used are - can we get these on a CSV file for reading/writing?
m = MimiBRICK.get_model()
run(m)

##==============================================================================
## DOECLIM+BRICK

#TODO - check what the parameter values used are - can we get these on a CSV file for reading/writing?
#function create_brick_doeclim(;rcp_scenario::String = "RCP85", start_year::Int=1850, end_year::Int=2300)
m = MimiBRICK_DOECLIM.create_brick_doeclim()
run(m)

##==============================================================================
## SNEASY+BRICK

#TODO - check what the parameter values used are - can we get these on a CSV file for reading/writing?
#function create_sneasy_brick(rcp_scenario::String; end_year::Int=2020)
m = create_sneasy_brick("RCP85")
run(m)

##==============================================================================
