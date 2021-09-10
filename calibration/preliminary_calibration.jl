##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
## This file carries out a preliminary Markov chain Monte Carlo calibration of BRICK.
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

# Model configuration
# --> Possible options: (1) "brick", (2) "doeclimbrick", (3) "sneasybrick"
model_config = "sneasybrick"

# Initial conditions from a previous file or from the prior distributions?
# --> If you want to use the midpoints of the prior ranges as the starting parameter estimates, and 5% of the prior range width as the step size for MCMC, set `start_from_priors = true`
start_from_priors = false
path_parameter_info = joinpath(@__DIR__, "..", "data", "calibration_data", "calibration_initial_values_"*model_config*".csv")
# --> If you want to read from a previous run, set these two file names/paths.
#     NOTE that if `start_from_priors = true`, these will NOT be used, even if they are set appropriately.
#     Also, the `path_initial_parameters` does not need to be distinct from the `path_parameter_info`.
#     `path_parameter_info` is just to get the prior ranges and names of the parameters, whereas `path_initial_parameters` will give the starting values for the model parameters.
if ~start_from_priors
    path_initial_parameters = joinpath(@__DIR__, "..", "data", "calibration_data", "calibration_initial_values_"*model_config*".csv")
    #     Also, the `path_initial_parameter_values` does not need to be distinct from the `path_initial_parameters`.
    #     `path_initial_parameters` is just
    path_initial_covariance = joinpath(@__DIR__, "..", "data", "calibration_data", "initial_proposal_covariance_matrix_"*model_config*".csv")
end

# Set years for model calibration.
calibration_start_year = 1850
calibration_end_year = 2017

# The length of the chain before burn-in and thinning
total_chain_length = 100_000
#total_chain_length = 100 # original was 100_000; this is for testing

# Burn-in length - How many samples from the beginning to immediately discard
# --> Not including as much burn-in as we would normally because the initial values are from the end of a 4-million iteration preliminary chain
burnin_length = 10_000

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


##------------------------------------------------------------------------------
## Initial set up
##------------------------------------------------------------------------------

# A folder with this name will be created to store all of the replication results.
results_folder_name = "my_"*model_config*"_results_"*string(Int(total_chain_length/1000))*"K_$(Dates.format(now(),"dd-mm-yyyy"))"

# Create output folder path for convenience and make path.
output = joinpath(@__DIR__, "..", "results", results_folder_name)
mkpath(output)

# Load calibration helper functions file.
include(joinpath("..", "calibration", "calibration_helper_functions.jl"))


##------------------------------------------------------------------------------
## Set up initial parameters and proposal covariance matrix
##------------------------------------------------------------------------------

# Either way, need the parameter names and distribution ranges
parameter_info = DataFrame(load(path_parameter_info, skiplines_begin=6))
num_parameters = nrow(parameter_info)
parnames = [Symbol(parameter_info.parameter[i]) for i in 1:num_parameters]
if start_from_priors
    # midpoint and 5% of the prior range width as the starting values and step sizes (step size ~ stdev, so squared to get variance)
    initial_parameters = 0.5*(parameter_info[:,"upper_bound"]+parameter_info[:,"lower_bound"])
    initial_covariance_matrix = diagm(0.05*(parameter_info[:,"upper_bound"]-parameter_info[:,"lower_bound"]))^2
else
    initial_parameters = DataFrame(load(path_initial_parameters, skiplines_begin=6)).starting_point
    initial_covariance_matrix = Array(Hermitian(Matrix(DataFrame(load(path_initial_covariance)))))
end


##------------------------------------------------------------------------------
## Load functions for running and calibrating the model configuration
## --> New configurations will need new drivers and posterior distribution calculation scripts, but can follow the format of the examples here
##------------------------------------------------------------------------------

# Load run historic model file.
include(joinpath("..", "calibration", "run_historic_models", "run_"*model_config*"_historic_climate.jl"))

# Load log-posterior script for this model configuration.
include(joinpath("..", "calibration", "create_log_posterior_"*model_config*".jl"))

# NOTE: the following two commands assume a naming convention for the `construct_run_[model_config]`
# and `construct_[model_config]_log_posterior` functions in the helper scripts included above.

# Create model run function used in log-posterior calculations.
@eval run_mymodel! = ($(Symbol("construct_run_$model_config")))(calibration_start_year, calibration_end_year)

# Create log-posterior function.
@eval log_posterior_mymodel = ($(Symbol("construct_$model_config"*"_log_posterior")))(run_mymodel!, model_start_year=1850, calibration_end_year=calibration_end_year, joint_antarctic_prior=false)

println("Begin baseline calibration of "*model_config*" model.\n")

# Carry out Bayesian calibration using robust adaptive metropolis MCMC algorithm.
Random.seed!(2021) # for reproducibility
#chain_raw, accept_rate, cov_matrix = RAM_sample(log_posterior_mymodel, initial_parameters, initial_covariance_matrix, Int(total_chain_length), opt_α=0.234)
chain_raw, accept_rate, cov_matrix, log_post = RAM_sample(log_posterior_mymodel, initial_parameters, initial_covariance_matrix, Int(total_chain_length), opt_α=0.234, output_log_probability_x=true)


##------------------------------------------------------------------------------
## Burn-in removal and check convergence via Gelman and Rubin potential scale
## reduction factor (PSRF)
##------------------------------------------------------------------------------

# Remove the burn-in period
chain_burned = chain_raw[(burnin_length+1):total_chain_length,:]

# Check convergence by computing Gelman and Rubin diagnostic for each parameter (potential scale reduction factor)
psrf = Array{Float64,1}(undef , num_parameters)
for p in 1:num_parameters
    chains = reshape(chain_burned[:,p], Int(size(chain_burned)[1]/num_walkers), num_walkers)
    chains = [chains[:,k] for k in 1:num_walkers]
    psrf[p] = potential_scale_reduction(chains...)
end

# Check if psrf < threshold_gr for each parameters
if all(x -> x < threshold_gr, psrf)
    println("All parameter chains have Gelman and Rubin PSRF < ",threshold_gr)
else
    println("WARNING: some parameter chains have Gelman and Rubin PSRF > ",threshold_gr)
    for p in 1:num_parameters
        println(parnames[p],"  ",round(psrf[p],digits=4))
    end
end

# this is just an example chunk of code for continuing a chain. to be deleted?
if false
    ## continuing a chain:
    # use the last iteration from the previous chain as the initial parameters.
    # use the last covariance matrix (that is output from RAM_sample) as the initial proposal covariance matrix.
    # run more iterations:
    num_new_samples = total_chain_length - burnin_length
    chain_raw_new, accept_rate_new, cov_matrix_new = RAM_sample(log_posterior_mymodel, chain_raw[total_chain_length,:], cov_matrix, Int(num_new_samples), opt_α=0.234)
    # reset to combine
    chain_raw = vcat(chain_raw, chain_raw_new)
    accept_rate = (total_chain_length*accept_rate + num_new_samples*accept_rate_new)/(total_chain_length+num_new_samples)
    cov_matrix = cov_matrix_new
end


##------------------------------------------------------------------------------
## Thinning - probably not done by default, just running many iterations and then subsampling
##------------------------------------------------------------------------------

if do_thinning
    # --> start the thinning lag out at the minimum `lags`
    # --> loop through the parameters and calculate the lag at which autocorrelation < threshold_acf
    # --> increase thinning lag whenever you encounter a higher needed lag
    thinlag = minimum(lags) # initialize
    for p in 1:num_parameters
        acf = autocor(chain_burned[:,p], lags)
        idx_low_enough = findall(x -> x <= threshold_acf, acf)
        if length(idx_low_enough) >= 1
            idx_low_enough = idx_low_enough[1]
        else
            longer_lags = lags[1]:(lags[2]-lags[1]):Int(0.5*(total_chain_length-burnin_length)) # increase the max lags
            acf = autocor(chain_burned[:,p], longer_lags)
            idx_low_enough = findall(x -> x <= threshold_acf, acf)[1]
            # TODO - the above will error out if there are still no lags at which ACF < threshold_acf
        end
        thinlag = maximum([thinlag, lags[idx_low_enough]])
    end
    # TODO - subsample the chain
    final_chain = DataFrame(chain_burned, :auto)
    rename!(final_chain, parnames)
else
    final_chain = DataFrame(chain_burned, :auto)
    rename!(final_chain, parnames)
end


##------------------------------------------------------------------------------
## Subsampling the final chains
##------------------------------------------------------------------------------

idx_subsample = sample(1:size(final_chain)[1], size_subsample, replace=false)
final_sample = final_chain[idx_subsample,:]


##------------------------------------------------------------------------------
## Save the results
##------------------------------------------------------------------------------

# Save calibrated parameter samples
println("Saving calibrated parameters for "*model_config*".\n")

save(joinpath(@__DIR__, output, "mcmc_acceptance_rate.csv"), DataFrame(acceptance_rate=accept_rate))
save(joinpath(@__DIR__, output, "proposal_covariance_matrix.csv"), DataFrame(cov_matrix, :auto))
save(joinpath(@__DIR__, output, "parameters_full_chain.csv"), final_chain)
save(joinpath(@__DIR__, output, "parameters_subsample.csv"), final_sample)

# Save initial conditions for future runs
path_new_initial_conditions = joinpath(@__DIR__, "..", "data", "calibration_data", "from_preliminary_chains")
mkpath(path_new_initial_conditions)
filename_new_initial_parameters = "calibration_initial_values_"*model_config*"_"*string(Int(total_chain_length/1000))*"K_$(Dates.format(now(),"dd-mm-yyyy")).csv"
new_initial_parameters = DataFrame(parameter_names = parnames, parameter_values = Vector(final_chain[size(final_chain)[1],:]))
save(joinpath(path_new_initial_conditions, filename_new_initial_parameters), new_initial_parameters)
filename_new_initial_covariance = "initial_proposal_covariance_matrix_"*model_config*"_"*string(Int(total_chain_length/1000))*"K_$(Dates.format(now(),"dd-mm-yyyy")).csv"
save(joinpath(path_new_initial_conditions, filename_new_initial_covariance), DataFrame(cov_matrix, :auto))


##------------------------------------------------------------------------------
## End
##------------------------------------------------------------------------------
