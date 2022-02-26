##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
## This file carries out a Markov chain Monte Carlo calibration of BRICK.
## This includes one of the following possible model configurations:
## (1) BRICK standalone (forced by input global mean surface temperatures and ocean heat uptake data)
## (2) DOECLIM+BRICK
## (3) SNEASY+BRICK
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------

# Activate the project for the paper and make sure all packages we need are installed.
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()


# Load required Julia packages.
using CSVFiles
using DataFrames
using Distributions
using KernelDensity
using LinearAlgebra
using Mimi
using NetCDF
using RobustAdaptiveMetropolisSampler
using MCMCDiagnostics
using Random
using StatsBase
using Dates

include(joinpath(@__DIR__, "calibration.jl"))
