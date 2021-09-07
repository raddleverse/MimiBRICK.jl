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
using Random
using StatsBase

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

# Do thinning?
# --> There are good reasons for thinning (accounting for autocorrelation in the Markov chain samples) and good reasons not to (e.g. Link and Eaton 2012; Maceachern and Berliner 1994)
# --> If `false`, the `threshold_acf` and `lags` below are not used.
do_thinning = false

# Threshold for the autocorrelation function (thinning; higher -> more permissive of autocorrelation in our samples)
threshold_acf = 0.05

# Lags - Which lags to check for computing the ACF? increments that aren't 1 are only used here for speed
lags = 10:10:50000


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
# --> These are the final values from a preliminary calibration 4 million iterations
# --> That preliminary calibration started from a SNEASY-BRICK calibration, but with the non-BRICK parameters removed
initial_parameters_brick = DataFrame(load(joinpath(@__DIR__, "..", "data", "calibration_data", "calibration_initial_values_brick.csv"), skiplines_begin=6))

# Load initial proposal covariance matrix (from previous calibrations) and format so it works with RAM sampler (need to account for rounding errors or Cholesky factorization fails).
# --> From the same preliminary calibration as the initial parameters above
initial_covariance_matrix_brick = Array(Hermitian(Matrix(DataFrame(load(joinpath(@__DIR__, "..", "data", "calibration_data", "initial_proposal_covariance_matrix_brick.csv"))))))

# Create `BRICK` function used in log-posterior calculations.
run_brick! = construct_run_brick(calibration_start_year, calibration_end_year)

# Create log-posterior function.
log_posterior_brick = construct_brick_log_posterior(run_brick!, model_start_year=1850, calibration_end_year=calibration_end_year, joint_antarctic_prior=false)

println("Begin baseline calibration of BRICK model.\n")

# Carry out Bayesian calibration using robust adaptive metropolis MCMC algorithm.
Random.seed!(2021) # for reproducibility
chain_brick, accept_rate_brick, cov_matrix_brick = RAM_sample(log_posterior_brick, initial_parameters_brick.starting_point, initial_covariance_matrix_brick, Int(total_chain_length), opt_α=0.234)

#-------------------------------------------------------------#
#-------------------------------------------------------------#
#------------------ Convergence checks and -------------------#
#------------------  (possibly) thinning   -------------------#
#-------------------------------------------------------------#
#-------------------------------------------------------------#

# Remove the burn-in period
chain_brick_burned = chain_brick[(burnin_length+1):total_chain_length,:]

# Check convergence by computing Gelman and Rubin diagnostic for each parameter (potential scale reduction factor)
num_parameters = size(chain_brick_burned)[2]
psrf = Array{Float64,1}(undef , num_parameters)
for p in 1:num_parameters
    chains = reshape(chain_brick_burned[:,p], Int(size(chain_brick_burned)[1]/num_walkers), num_walkers)
    chains = [chains[:,k] for k in 1:num_walkers]
    psrf[p] = potential_scale_reduction(chains...)
end

# Check if psrf < threshold_gr for each parameters
if all(x -> x < threshold_gr, psrf)
    println("All parameter chains have Gelman and Rubin PSRF < ",threshold_gr)
else
    println("WARNING: some parameter chains have Gelman and Rubin PSRF > ",threshold_gr)
    for p in 1:num_parameters
        println(initial_parameters_brick.parameter[p],"  ",round(psrf[p],digits=4))
    end
end

if false
    ## continuing a chain:
    # set the initial parameter sample
    initial_parameters_brick[:,"starting_point"] = chain_brick[total_chain_length,:]
    # run more iterations
    num_new_samples = total_chain_length - burnin_length
    chain_brick2, accept_rate_brick2, cov_matrix_brick2 = RAM_sample(log_posterior_brick, initial_parameters_brick.starting_point, cov_matrix_brick, Int(num_new_samples), opt_α=0.234)
    full_chain_brick = vcat(chain_brick, chain_brick2)
end


# Thinning
if do_thinning
    # --> start the thinning lag out at the minimum `lags`
    # --> loop through the parameters and calculate the lag at which autocorrelation < threshold_acf
    # --> increase thinning lag whenever you encounter a higher needed lag
    thinlag = minimum(lags) # initialize
    for p in 1:num_parameters
        acf = autocor(chain_brick_burned[:,p], lags)
        idx_low_enough = findall(x -> x <= threshold_acf, acf)
        if length(idx_low_enough) >= 1
            idx_low_enough = idx_low_enough[1]
        else
            longer_lags = lags[1]:(lags[2]-lags[1]):Int(0.5*(total_chain_length-burnin_length)) # increase the max lags
            acf = autocor(chain_brick_burned[:,p], longer_lags)
            idx_low_enough = findall(x -> x <= threshold_acf, acf)[1]
            # TODO - the above will error out if there are still no lags at which ACF < threshold_acf
        end
        thinlag = maximum([thinlag, lags[idx_low_enough]])
    end
    # TODO - subsample the chain
    final_chain_brick = DataFrame(chain_brick_burned, :auto)
    rename!(final_chain_brick, [Symbol(initial_parameters_brick.parameter[i]) for i in 1:num_parameters])
else
    final_chain_brick = DataFrame(chain_brick_burned, :auto)
    rename!(final_chain_brick, [Symbol(initial_parameters_brick.parameter[i]) for i in 1:num_parameters])
end

# Subsample
idx_subsample = sample(1:size(final_chain_brick)[1], size_subsample, replace=false)
final_sample_brick = final_chain_brick[idx_subsample,:]

#--------------------------------------------------#
#------------ Save Calibration Results ------------#
#--------------------------------------------------#

# Save calibrated parameter samples
println("Saving calibrated parameters for BRICK.\n")

# BRICK model calibration.
save(joinpath(@__DIR__, output, "mcmc_acceptance_rate.csv"), DataFrame(brick_acceptance=accept_rate_brick))
save(joinpath(@__DIR__, output, "proposal_covariance_matrix.csv"), DataFrame(cov_matrix_brick, :auto))
save(joinpath(@__DIR__, output, "parameters_full_chain.csv"), final_chain_brick)
save(joinpath(@__DIR__, output, "parameters_subsample.csv"), final_sample_brick)
