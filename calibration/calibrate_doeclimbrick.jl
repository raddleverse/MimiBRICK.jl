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
using MCMCDiagnostics
using Random
using StatsBase

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

# The length of the chain before burn-in and thinning
total_chain_length = 11_000_000
#total_chain_length = 100 # original was 100_000; this is for testing

# Burn-in length - How many samples from the beginning to immediately discard
# --> Not including as much burn-in as we would normally because the initial values are from the end of a 4-million iteration preliminary chain
burnin_length = 1_000_000

# Threshold for Gelman and Rubin potential scale reduction factor (burn-in/convergence)
# --> 1.1 or 1.05 are standard practice. Further from 1 is
threshold_gr = 1.1

# Number of sub-chains that the single larger chain is divided into for convergence check
# --> Convergence check will chop off the burn-in period from the beginning, then divide the remaining chain into `num_walkers` equal-length chains. Then
# --> the GR diagnostic (PSRF) will check how the within-chain variance compares to the between-chain variance. Converged if they're similar.
# --> Note that fewer walkers will also be more permissive, since (all other things being equal) variance is higher with fewer samples (# walkers)
num_walkers = 2

# Create a subsample of posterior parameters?
size_subsample = 10_000


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

# Load inital parameter values for BRICK model.
# --> These are the final values from a preliminary calibration 4 million iterations
# --> That preliminary calibration started from a SNEASY-BRICK calibration, but with the non-DOECLIM/BRICK parameters removed
initial_parameters_doeclimbrick = DataFrame(load(joinpath(@__DIR__, "..", "data", "calibration_data", "calibration_initial_values_doeclim_brick.csv"), skiplines_begin=6))
num_parameters = nrow(initial_parameters_doeclimbrick)
parnames = [Symbol(initial_parameters_doeclimbrick.parameter[i]) for i in 1:num_parameters]

# Load initial proposal covariance matrix (from previous calibrations) and format so it works with RAM sampler (need to account for rounding errors or Cholesky factorization fails).
# --> From the same preliminary calibration as the initial parameters above
initial_covariance_matrix_doeclimbrick = Array(Hermitian(Matrix(DataFrame(load(joinpath(@__DIR__, "..", "data", "calibration_data", "initial_proposal_covariance_matrix_doeclimbrick.csv"))))))

# Create `DOECLIM+BRICK` function used in log-posterior calculations.
run_doeclimbrick! = construct_run_doeclimbrick(calibration_start_year, calibration_end_year)

# Create log-posterior function.
log_posterior_doeclim_brick = construct_doeclimbrick_log_posterior(run_doeclimbrick!, model_start_year=1850, calibration_end_year=calibration_end_year, joint_antarctic_prior=false)

println("Begin baseline calibration of DOECLIM+BRICK model.\n")

# Carry out Bayesian calibration using robust adaptive metropolis MCMC algorithm.
Random.seed!(2021) # for reproducibility
chain_doeclimbrick, accept_rate_doeclimbrick, cov_matrix_doeclimbrick = RAM_sample(log_posterior_doeclim_brick, initial_parameters_doeclimbrick.starting_point, initial_covariance_matrix_doeclimbrick, Int(total_chain_length), opt_α=0.234)

#-------------------------------------------------------------#
#-------------------------------------------------------------#
#------------------ Convergence checks and -------------------#
#------------------  (possibly) thinning   -------------------#
#-------------------------------------------------------------#
#-------------------------------------------------------------#

# Remove the burn-in period
chain_doeclimbrick_burned = chain_doeclimbrick[(burnin_length+1):total_chain_length,:]

# Check convergence by computing Gelman and Rubin diagnostic for each parameter (potential scale reduction factor)
psrf = Array{Float64,1}(undef , num_parameters)
for p in 1:num_parameters
    chains = reshape(chain_doeclimbrick_burned[:,p], Int(size(chain_doeclimbrick_burned)[1]/num_walkers), num_walkers)
    chains = [chains[:,k] for k in 1:num_walkers]
    psrf[p] = potential_scale_reduction(chains...)
end

# Check if psrf < threshold_gr for each parameters
if all(x -> x < threshold_gr, psrf)
    println("All parameter chains have Gelman and Rubin PSRF < ",threshold_gr)
else
    println("WARNING: some parameter chains have Gelman and Rubin PSRF > ",threshold_gr)
    for p in 1:num_parameters
        println(initial_parameters_doeclimbrick.parameter[p],"  ",round(psrf[p],digits=4))
    end
end

# Thinning? By default, not thinning for BRICK. Can add here though if desired.
final_chain_doeclimbrick = DataFrame(chain_doeclimbrick_burned, :auto)
rename!(final_chain_doeclimbrick, [Symbol(initial_parameters_doeclimbrick.parameter[i]) for i in 1:num_parameters])

# Subsample
idx_subsample = sample(1:size(final_chain_doeclimbrick)[1], size_subsample, replace=false)
final_sample_doeclimbrick = final_chain_doeclimbrick[idx_subsample,:]

#--------------------------------------------------#
#------------ Save Calibration Results ------------#
#--------------------------------------------------#

# Save calibrated parameter samples
println("Saving calibrated parameters for DOECLIM+BRICK.\n")

# DOECLIM-BRICK model calibration.
save(joinpath(@__DIR__, output, "mcmc_acceptance_rate.csv"), DataFrame(doeclimbrick_acceptance=accept_rate_doeclimbrick))
save(joinpath(@__DIR__, output, "proposal_covariance_matrix.csv"), DataFrame(cov_matrix_doeclimbrick, :auto))
save(joinpath(@__DIR__, output, "parameters_full_chain.csv"), final_chain_doeclimbrick)
save(joinpath(@__DIR__, output, "parameters_subsample.csv"), final_sample_doeclimbrick)
