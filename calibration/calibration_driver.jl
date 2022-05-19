using MimiBRICK
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

##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
## This file carries out a Markov chain Monte Carlo calibration of BRICK.
## This includes one of the following possible model configurations:
## (1) BRICK standalone (forced by input global mean surface temperatures and ocean heat uptake data)
## (2) DOECLIM+BRICK
## (3) SNEASY+BRICK
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------

calibration_start_year = 1850
calibration_end_year   = 2017
total_chain_length     = 20_000_000
size_subsample         = 10_000
threshold_gr           = 1.1

## Create the log-posterior functions
log_posterior_brick = MimiBRICK.construct_brick_log_posterior(MimiBRICK.construct_run_brick(calibration_start_year, calibration_end_year), model_start_year=calibration_start_year, calibration_end_year=calibration_end_year, joint_antarctic_prior=false)
log_posterior_doeclimbrick = MimiBRICK.construct_doeclimbrick_log_posterior(MimiBRICK.construct_run_doeclimbrick(calibration_start_year, calibration_end_year), model_start_year=calibration_start_year, calibration_end_year=calibration_end_year, joint_antarctic_prior=false)
log_posterior_sneasybrick = MimiBRICK.construct_sneasybrick_log_posterior(MimiBRICK.construct_run_sneasybrick(calibration_start_year, calibration_end_year), model_start_year=calibration_start_year, calibration_end_year=calibration_end_year, joint_antarctic_prior=false)

## Do the actual calibrations

# BRICK calibration
x = MimiBRICK.run_calibration(log_posterior_brick; model_config="brick", calibration_start_year=calibration_start_year,
                    calibration_end_year=calibration_end_year, total_chain_length=total_chain_length,
                    burnin_length=1_000_000, threshold_gr=threshold_gr, num_walkers=2,
                    size_subsample=size_subsample, start_from_priors=false)

# DOECLIM-BRICK calibration
x = MimiBRICK.run_calibration(log_posterior_doeclimbrick; model_config="doeclimbrick", calibration_start_year=1850, calibration_end_year=2017,
                    total_chain_length=20_000_000, burnin_length=7_000_000, threshold_gr=1.1, num_walkers=2,
                    size_subsample=10_000, start_from_priors=false)

# SNEASY-BRICK calibration
x = MimiBRICK.run_calibration(model_config="sneasybrick", calibration_start_year=1850, calibration_end_year=2017,
                    total_chain_length=20_000_000, burnin_length=1_000_000, threshold_gr=1.1, num_walkers=2,
                    size_subsample=10_000, start_from_priors=false)
