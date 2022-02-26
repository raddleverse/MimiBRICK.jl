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

# BRICK calibration
x = run_calibration(model_config="brick", calibration_start_year=1850, calibration_end_year=2017,
                    total_chain_length=20_000_000, burnin_length=1_000_000, threshold_gr=1.1, num_walkers=2,
                    size_subsample=10_000, start_from_priors=false)

# DOECLIM-BRICK calibration
x = run_calibration(model_config="doeclimbrick", calibration_start_year=1850, calibration_end_year=2017,
                    total_chain_length=20_000_000, burnin_length=7_000_000, threshold_gr=1.1, num_walkers=2,
                    size_subsample=10_000, start_from_priors=false)

# SNEASY-BRICK calibration
x = run_calibration(model_config="sneasybrick", calibration_start_year=1850, calibration_end_year=2017,
                    total_chain_length=20_000_000, burnin_length=1_000_000, threshold_gr=1.1, num_walkers=2,
                    size_subsample=10_000, start_from_priors=false)
