##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
## This function carries out a Markov chain Monte Carlo calibration of BRICK.
## This includes one of the following possible model configurations:
## (1) BRICK standalone (forced by input global mean surface temperatures and ocean heat uptake data)
## (2) DOECLIM+BRICK
## (3) SNEASY+BRICK
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------

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


function run_calibration(log_posterior_mymodel; model_config="brick", calibration_start_year=1850, calibration_end_year=2005,
                         total_chain_length=1000, burnin_length=0, threshold_gr=1.1, num_walkers=2,
                         size_subsample=1000, start_from_priors=false)

    # File name and path to obtain the parameter names (*in order*)
    path_parameter_info = joinpath(@__DIR__, "..", "data", "calibration_data", "calibration_initial_values_"*model_config*".csv")
    # If you want to read from a previous run, set these two file names/paths.
    # NOTE that if `start_from_priors = true`, these will NOT be used, even if they are set appropriately.
    # Also, the `path_initial_parameters` does not need to be distinct from the `path_parameter_info`.
    # `path_parameter_info` is just to get the names of the parameters, whereas `path_initial_parameters` will provide the starting values for the model parameters as well.
    if ~start_from_priors
        path_initial_parameters = joinpath(@__DIR__, "..", "data", "calibration_data", "calibration_initial_values_"*model_config*".csv")
        path_initial_covariance = joinpath(@__DIR__, "..", "data", "calibration_data", "initial_proposal_covariance_matrix_"*model_config*".csv")
    end

    ##------------------------------------------------------------------------------
    ## Initial set up
    ##------------------------------------------------------------------------------

    # A folder with this name will be created to store all of the replication results.
    if total_chain_length < 1_000 # hundreds of iterations
        chain_len_str = string(Int(total_chain_length))
    elseif total_chain_length < 1_000_000 # thousands of iterations
        chain_len_str = string(Int(total_chain_length/1000))*"K"
    elseif total_chain_length < 1_000_000_000 # millions of iterations
        chain_len_str = string(Int(total_chain_length/1000000))*"M"
    else # billions of iterations (or more, but probably not)
        chain_len_str = string(Int(total_chain_length/1000000000))*"B"
    end
    results_folder_name = "my_"*model_config*"_results_"*chain_len_str*"_$(Dates.format(now(),"dd-mm-yyyy"))"

    # Create output folder path for convenience and make path.
    output = joinpath(@__DIR__, "..", "results", results_folder_name)
    mkpath(output)

    ##------------------------------------------------------------------------------
    ## Set up initial parameters and proposal covariance matrix
    ##------------------------------------------------------------------------------

    # Either way, need the parameter names and distribution ranges
    parameter_info = DataFrame(load(path_parameter_info))
    num_parameters = nrow(parameter_info)
    parnames = [Symbol(parameter_info.parameter_names[i]) for i in 1:num_parameters]
    if start_from_priors
        error("Starting from priors not supported yet.")
        # midpoint and 5% of the prior range width as the starting values and step sizes (step size ~ stdev, so squared to get variance)
        initial_parameters = 0.5*(parameter_info[:,"upper_bound"]+parameter_info[:,"lower_bound"])
        initial_covariance_matrix = diagm(0.05*(parameter_info[:,"upper_bound"]-parameter_info[:,"lower_bound"]))^2
    else
        initial_parameters = DataFrame(load(path_initial_parameters)).parameter_values
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

    # NOTE: the following two commands for each model config assume a naming convention
    # for the `construct_run_[model_config]` and `construct_[model_config]_log_posterior`
    # functions in the helper scripts included above. Using this instead of the
    # @eval and Symbols so this can be run as a function instead of a script.
    if false
    if model_config=="brick"
        run_mymodel! = construct_run_brick(calibration_start_year, calibration_end_year)
        log_posterior_mymodel = construct_brick_log_posterior(run_mymodel!, model_start_year=calibration_start_year, calibration_end_year=calibration_end_year, joint_antarctic_prior=false)
    elseif model_config=="doeclimbrick"
        run_mymodel! = construct_run_doeclimbrick(calibration_start_year, calibration_end_year)
        log_posterior_mymodel = construct_doeclimbrick_log_posterior(run_mymodel!, model_start_year=calibration_start_year, calibration_end_year=calibration_end_year, joint_antarctic_prior=false)
    elseif model_config=="sneasybrick"
        run_mymodel! = construct_run_sneasybrick(calibration_start_year, calibration_end_year)
        log_posterior_mymodel = construct_sneasybrick_log_posterior(run_mymodel!, model_start_year=calibration_start_year, calibration_end_year=calibration_end_year, joint_antarctic_prior=false)
    end
end

    println("Begin baseline calibration of "*model_config*" model.\n")

    # Carry out Bayesian calibration using robust adaptive metropolis MCMC algorithm.
    Random.seed!(2021) # for reproducibility
    @time chain_raw, accept_rate, cov_matrix, log_post = RAM_sample(log_posterior_mymodel, initial_parameters, initial_covariance_matrix, Int(total_chain_length), opt_α=0.234, output_log_probability_x=true)

    ##------------------------------------------------------------------------------
    ## Burn-in removal and check convergence via Gelman and Rubin potential scale
    ## reduction factor (PSRF)
    ##------------------------------------------------------------------------------

    # Remove the burn-in period
    chain_burned = chain_raw[(burnin_length+1):total_chain_length,:]
    log_post_burned = log_post[(burnin_length+1):total_chain_length]

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


    ##------------------------------------------------------------------------------
    ## Subsampling the final chains
    ##------------------------------------------------------------------------------

    Random.seed!(2022) # for reproducibility
    idx_subsample = sample(1:size(chain_burned)[1], size_subsample, replace=false)
    final_sample = chain_burned[idx_subsample,:]
    log_post_final_sample = log_post_burned[idx_subsample]

    ##------------------------------------------------------------------------------
    ## Save the results
    ##------------------------------------------------------------------------------

    # Save calibrated parameter samples
    println("Saving calibrated parameters for "*model_config*".\n")

    save(joinpath(@__DIR__, output, "mcmc_log_post.csv"), DataFrame(log_post=log_post))
    save(joinpath(@__DIR__, output, "mcmc_acceptance_rate.csv"), DataFrame(acceptance_rate=accept_rate))
    save(joinpath(@__DIR__, output, "proposal_covariance_matrix.csv"), DataFrame(cov_matrix, :auto))
    save(joinpath(@__DIR__, output, "parameters_full_chain.csv"), DataFrame(chain_raw,parnames))
    save(joinpath(@__DIR__, output, "parameters_subsample.csv"), DataFrame(final_sample,parnames))
    save(joinpath(@__DIR__, output, "log_post_subsample.csv"), DataFrame(log_post=log_post_final_sample))

    # Save initial conditions for future runs
    path_new_initial_conditions = joinpath(@__DIR__, "..", "data", "calibration_data", "from_calibration_chains")
    mkpath(path_new_initial_conditions)
    filename_new_initial_parameters = "calibration_initial_values_"*model_config*"_"*chain_len_str*"_$(Dates.format(now(),"dd-mm-yyyy")).csv"
    new_initial_parameters = DataFrame(parameter_names = parnames, parameter_values = Vector(chain_burned[size(chain_burned)[1],:]))
    save(joinpath(path_new_initial_conditions, filename_new_initial_parameters), new_initial_parameters)
    filename_new_initial_covariance = "initial_proposal_covariance_matrix_"*model_config*"_"*chain_len_str*"_$(Dates.format(now(),"dd-mm-yyyy")).csv"
    save(joinpath(path_new_initial_conditions, filename_new_initial_covariance), DataFrame(cov_matrix, :auto))

    return (DataFrame(chain_raw,parnames), accept_rate, cov_matrix, log_post, DataFrame(final_sample,parnames), log_post_final_sample)
end

##------------------------------------------------------------------------------
## End
##------------------------------------------------------------------------------
