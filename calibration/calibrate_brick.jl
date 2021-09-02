# #-------------------------------------------------------------------------------------------------------
# #-------------------------------------------------------------------------------------------------------
# # This file carries out a Markov chain Monte Carlo calibration of BRICK.
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
using MCMCDiagnostics


# A folder with this name will be created to store all of the replication results.
results_folder_name = "my_brick_results"

# Create output folder path for convenience and make path.
output = joinpath(@__DIR__, "..", "results", results_folder_name)
mkpath(output)

# Load calibration helper functions file.
include(joinpath("..", "calibration", "calibration_helper_functions.jl"))

# Set years for model calibration.
calibration_start_year = 1850
calibration_end_year = 2017

# The length of the chain before burn-in and thinning
total_chain_length = 2_000_000
#total_chain_length = 100 # original was 100_000; this is for testing

# Burn-in length - How many samples from the beginning to immediately discard
burnin_length = 1_000_000

# Number of individual hchain_brick_burned = chain_brick[(burnin_length+1)ypothetical chains to break the big on up into for testing convergence ("walkers"):total_chain_length,:]

chain_brick_burned = chain_brick[(burnin_length+1):total_chain_length,:]


#-------------------------------------------------------------#
#-------------------------------------------------------------#
#---------------- BRICK Baseline Calibration -----------------#
#-------------------------------------------------------------#
#-------------------------------------------------------------#

# NOTE** This version uses the kernel density estimated marginal priors for the Antarctic ice sheet based on a calibration to paleo data.

# Load run historic model file.
include(joinpath("..", "calibration", "run_historic_models", "run_brick_historic_climate.jl"))

# Load log-posterior script for BRICK model.
include(joinpath("..", "calibration", "create_log_posterior_brick.jl"))

# Load inital parameter values for BRICK model.
# These are from a BRICK calibration, but with the non-BRICK parameters removed
initial_parameters_brick = DataFrame(load(joinpath(@__DIR__, "..", "data", "calibration_data", "calibration_initial_values_brick.csv"), skiplines_begin=6))

# Load initial proposal covariance matrix (from previous calibrations) and format so it works with RAM sampler (need to account for rounding errors or Cholesky factorization fails).
initial_covariance_matrix_brick = Array(Hermitian(Matrix(DataFrame(load(joinpath(@__DIR__, "..", "data", "calibration_data", "initial_proposal_covariance_matrix_brick.csv"))))))

# Create `BRICK` function used in log-posterior calculations.
run_brick! = construct_run_brick(calibration_start_year, calibration_end_year)

# Create log-posterior function.
log_posterior_brick = construct_brick_log_posterior(run_brick!, model_start_year=1850, calibration_end_year=calibration_end_year, joint_antarctic_prior=false)

println("Begin baseline calibration of BRICK model.\n")

# Carry out Bayesian calibration using robust adaptive metropolis MCMC algorithm.
chain_brick, accept_rate_brick, cov_matrix_brick = RAM_sample(log_posterior_brick, initial_parameters_brick.starting_point, initial_covariance_matrix_brick, Int(total_chain_length), opt_α=0.234)

# No processing done here; doing that outside the calibration

# or maybe processing done here...??

# burn in?
chain_brick_burned = chain_brick[(burnin_length+1):total_chain_length,:]

# compute Gelman and Rubin diagnostic for each parameter
# (potential scale reduction factor)
num_parameters = size(chain_brick_burned)[2]
psrf = Array{Float64,1}(undef , num_parameters)
for p in 1:num_parameters
    chains = reshape(chain_brick_burned[:,p], Int(size(chain_brick_burned)[1]/num_walkers), num_walkers)
    chains = [chains[:,k] for k in 1:num_walkers]
    psrf[p] = potential_scale_reduction(chains...)
end

## continuing a chain:
# set the initial parameter sample
initial_parameters_brick[:,"starting_point"] = chain_brick[total_chain_length,:]
# run more iterations
num_new_samples = 2_000_000
chain_brick2, accept_rate_brick2, cov_matrix_brick2 = RAM_sample(log_posterior_brick, initial_parameters_brick.starting_point, cov_matrix_brick, Int(num_new_samples), opt_α=0.234)
full_chain_brick = vcat(chain_brick, chain_brick2)

#--------------------------------------------------#
#------------ Save Calibration Results ------------#
#--------------------------------------------------#

# Save calibrated parameter samples
println("Saving calibrated parameters for BRICK.\n")

# BRICK model calibration.
save(joinpath(@__DIR__, output, "mcmc_acceptance_rate.csv"), DataFrame(brick_acceptance=accept_rate_brick))
save(joinpath(@__DIR__, output, "proposal_covariance_matrix.csv"), DataFrame(cov_matrix_brick, :auto))
save(joinpath(@__DIR__, output, "parameters_full_chain.csv"), DataFrame(chain_brick, [Symbol(initial_parameters_brick.parameter[i]) for i in 1:length(mean_brick)]))
