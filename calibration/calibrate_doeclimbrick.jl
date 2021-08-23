# #-------------------------------------------------------------------------------------------------------
# #-------------------------------------------------------------------------------------------------------
# # This file carries out a Markov chain Monte Carlo calibration of DOECLIM+BRICK.
# #-------------------------------------------------------------------------------------------------------
# #-------------------------------------------------------------------------------------------------------

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


# A folder with this name will be created to store all of the replication results.
results_folder_name = "my_doeclim+brick_results"

# Create output folder path for convenience and make path.
output = joinpath(@__DIR__, "..", "results", results_folder_name)
mkpath(output)

# Load calibration helper functions file.
include(joinpath("..", "calibration", "calibration_helper_functions.jl"))

# Set years for model calibration.
calibration_start_year = 1850
calibration_end_year = 2017

# The length of the final chain (i.e. number of samples from joint posterior pdf after discarding burn-in period values).
#final_chain_length = 100_000
final_chain_length = 100 # original was 100_000; this is for testing

# Length of burn-in period (i.e. number of initial MCMC samples to discard).
#burn_in_length = 1_000
burn_in_length = 10 # original was 1_000; this is for testing


#-------------------------------------------------------------#
#-------------------------------------------------------------#
#------------ DOECLIM + BRICK Baseline Calibration -----------#
#-------------------------------------------------------------#
#-------------------------------------------------------------#

# NOTE** This version uses the kernel density estimated marginal priors for the Antarctic ice sheet based on a calibration to paleo data.

# Load run historic model file.
include(joinpath("..", "calibration", "run_historic_models", "run_doeclim_brick_historic_climate.jl"))

# Load log-posterior script for DOECLIM+BRICK model.
include(joinpath("..", "calibration", "create_log_posterior_doeclim_brick.jl"))

# Load inital parameter values for DOECLIM+BRICK model.
# These are from a DOECLIM-BRICK calibration, but with the non-DOECLIM parameters removed
initial_parameters_doeclimbrick = DataFrame(load(joinpath(@__DIR__, "..", "data", "calibration_data", "calibration_initial_values_doeclim_brick.csv"), skiplines_begin=6))

# Load initial proposal covariance matrix (from previous calibrations) and format so it works with RAM sampler (need to account for rounding errors or Cholesky factorization fails).
initial_covariance_matrix_doeclimbrick = Array(Hermitian(Matrix(DataFrame(load(joinpath(@__DIR__, "..", "data", "calibration_data", "initial_proposal_covariance_matrix_doeclimbrick.csv"))))))

# Create `DOECLIM+BRICK` function used in log-posterior calculations.
run_doeclimbrick! = construct_run_doeclimbrick(calibration_start_year, calibration_end_year)

# Create log-posterior function.
log_posterior_doeclim_brick = construct_doeclimbrick_log_posterior(run_doeclimbrick!, model_start_year=1850, calibration_end_year=calibration_end_year, joint_antarctic_prior=false)

println("Begin baseline calibration of DOECLIM+BRICK model.\n")

# Carry out Bayesian calibration using robust adaptive metropolis MCMC algorithm.
chain_doeclimbrick, accept_rate_doeclimbrick, cov_matrix_doeclimbrick = RAM_sample(log_posterior_doeclim_brick, initial_parameters_doeclimbrick.starting_point, initial_covariance_matrix_doeclimbrick, Int(final_chain_length + burn_in_length), opt_α=0.234)

# Discard burn-in values.
burned_chain_doeclimbrick = chain_doeclimbrick[Int(burn_in_length+1):end, :]

# Calculate mean posterior parameter values.
mean_doeclimbrick = vec(mean(burned_chain_doeclimbrick, dims=1))

# Calculate posterior correlations between parameters and set column names.
correlations_doeclimbrick = DataFrame(cor(burned_chain_doeclimbrick), :auto)
rename!(correlations_doeclimbrick, [Symbol(initial_parameters_doeclimbrick.parameter[i]) for i in 1:length(mean_doeclimbrick)])

# Create equally-spaced indices to thin chains down to 10,000 and 100,000 samples.
thin_indices_100k = trunc.(Int64, collect(range(1, stop=final_chain_length, length=100_000)))
thin_indices_10k  = trunc.(Int64, collect(range(1, stop=final_chain_length, length=10_000)))

# Create thinned chains (after burn-in period) with 10,000 and 100,000 samples and assign parameter names to each column.
thin100k_chain_doeclimbrick = DataFrame(burned_chain_doeclimbrick[thin_indices_100k, :], :auto)
thin10k_chain_doeclimbrick  = DataFrame(burned_chain_doeclimbrick[thin_indices_10k, :], :auto)

rename!(thin100k_chain_doeclimbrick, [Symbol(initial_parameters_doeclimbrick.parameter[i]) for i in 1:length(mean_doeclimbrick)])
rename!(thin10k_chain_doeclimbrick,  [Symbol(initial_parameters_doeclimbrick.parameter[i]) for i in 1:length(mean_doeclimbrick)])

#--------------------------------------------------#
#------------ Save Calibration Results ------------#
#--------------------------------------------------#

# Save calibrated parameter samples
println("Saving calibrated parameters for DOECLIM+BRICK.\n")

# DOECLIM-BRICK model calibration.
save(joinpath(@__DIR__, output, "mcmc_acceptance_rate.csv"), DataFrame(doeclimbrick_acceptance=accept_rate_doeclimbrick))
save(joinpath(@__DIR__, output, "proposal_covariance_matrix.csv"), DataFrame(cov_matrix_doeclimbrick, :auto))
save(joinpath(@__DIR__, output, "mean_parameters.csv"), DataFrame(parameter = initial_parameters_doeclimbrick.parameter[1:length(mean_doeclimbrick)], doeclimbrick_mean=mean_doeclimbrick))
save(joinpath(@__DIR__, output, "parameters_10k.csv"), thin10k_chain_doeclimbrick)
save(joinpath(@__DIR__, output, "parameters_100k.csv"), thin100k_chain_doeclimbrick)
save(joinpath(@__DIR__, output, "posterior_correlations.csv"), correlations_doeclimbrick)
