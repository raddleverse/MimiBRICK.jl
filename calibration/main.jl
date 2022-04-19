# #-------------------------------------------------------------------------------------------------------
# #-------------------------------------------------------------------------------------------------------
# # This file carries out all of the runs to replicate the results from "BRICK-SCC Paper".
# #-------------------------------------------------------------------------------------------------------
# #-------------------------------------------------------------------------------------------------------

# Load required Julia packages.
using CSVFiles
using DataFrames
using Distributions
using Interpolations
using KernelDensity
using LinearAlgebra
using Mimi
using NetCDF
using RobustAdaptiveMetropolisSampler

# A folder with this name will be created to store all of the replication results.
#results_folder_name = "my_results"

#####################################################################################################
#####################################################################################################
#                           Run Social Cost of Carbon Replication Code                             #
#####################################################################################################
#####################################################################################################

# Load generic helper functions file.
include("helper_functions.jl")

# Create folder structure to save results.
#build_result_folders(results_folder_name)

# Create output folder path for convenience.
#output = joinpath("..", "results", results_folder_name)

#######################################################################
#---------------------------------------------------------------------#
#                    Calibrate Climate Models                         #
#---------------------------------------------------------------------#
#######################################################################

println("BEGIN CALIBRATING MODELS.\n")

# Load calibration helper functions file common to all models.
include(joinpath("..", "calibration", "calibration_helper_functions.jl"))

# Set final year for model calibration.
calibration_end_year = 2017

# The length of the final chain (i.e. number of samples from joint posterior pdf after discarding burn-in period values).
final_chain_length = 2_500_000

# Length of burn-in period (i.e. number of initial MCMC samples to discard).
burn_in_length = final_chain_length * 0.2

# Select number of total samples (final samples + burn-in).
n_mcmc_samples = Int(final_chain_length + burn_in_length)

# Create equally-spaced indices to thin chains down to 10,000 and 100,000 samples.
thin_indices_100k = trunc.(Int64, collect(range(1, stop=final_chain_length, length=100_000)))
thin_indices_10k  = trunc.(Int64, collect(range(1, stop=final_chain_length, length=10_000)))


#---------------------------------------#
#---------------------------------------#
#------- PAGE Climate (with ECS) -------#
#---------------------------------------#
#---------------------------------------#

# Load run historic model file.
include(joinpath(@__DIR__, "..", "calibration", "run_historic_models", "run_page_historic_climate_ecs.jl"))

# Load log-posterior file for PAGE climate model.
include(joinpath(@__DIR__, "..", "calibration", "create_log_posterior_page_climate.jl"))

# Load inital parameter values for the PAGE climate model (with the equilibrium climate sensitivity set as a parameter).
initial_parameters = DataFrame(load(joinpath(@__DIR__, "..", "data", "calibration_data", "calibration_initial_values_page_climate.csv"), skiplines_begin=6))

# Create `run_page_climate_ecs` function used in log-posterior calculations.
run_page_climate! = construct_run_page_climate_ecs(calibration_end_year)

#----------------------------------#
#------ Baseline Calibration ------#
#----------------------------------#

println("Begin calibration of PAGE climate model with ECS parameter.\n")

# Create log-posterior function.
log_posterior_page_climate = construct_page_climate_log_posterior(run_page_climate!, model_start_year=1850, calibration_end_year=calibration_end_year)

# Carry out Bayesian calibration using robust adaptive metropolis MCMC algorithm.
chain_page_climate, accept_rate_page_climate, cov_matrix_page_climate = RAM_sample(log_posterior_page_climate, initial_parameters.starting_point, Diagonal(initial_parameters.mcmc_step_size), n_mcmc_samples, opt_α=0.234)

# Discard burn-in values.
burned_chain_page_climate = chain_page_climate[Int(burn_in_length+1):end, :]

# Calculate mean posterior parameter values.
mean_page_climate = vec(mean(burned_chain_page_climate, dims=1))

# Create thinned chains (after burn-in period) with 10,000 and 100,000 samples and assign parameter names to each column.
thin100k_chain_page_climate = DataFrame(burned_chain_page_climate[thin_indices_100k, :])
thin10k_chain_page_climate  = DataFrame(burned_chain_page_climate[thin_indices_10k, :])

names!(thin100k_chain_page_climate, [Symbol(initial_parameters.parameter[i]) for i in 1:length(mean_page_climate)])
names!(thin10k_chain_page_climate,  [Symbol(initial_parameters.parameter[i]) for i in 1:length(mean_page_climate)])

# !!! TO DO, SET UP OUTPUT FOLDER DIRECTORY FOR SAVING RESULTS !!!
output = joinpath("results", "page_climate_ecs_gmsl_params_outside_priorrange_may24")
mkdir(joinpath(output, "calibrated_parameters"))
save(joinpath(@__DIR__, output, "calibrated_parameters", "mcmc_acceptance_rate.csv"), DataFrame(page_acceptance=accept_rate_page_climate))
save(joinpath(@__DIR__, output, "calibrated_parameters", "mean_parameters.csv"), DataFrame(parameter = initial_parameters.parameter[1:length(mean_page_climate)], mean=mean_page_climate))
save(joinpath(@__DIR__, output, "calibrated_parameters", "parameters_10k.csv"), DataFrame(thin10k_chain_page_climate))
save(joinpath(@__DIR__, output, "calibrated_parameters", "parameters_100k.csv"), DataFrame(thin100k_chain_page_climate))



#-----------------------------------------#
#-----------------------------------------#
#------- DICE Climate (annualized) -------#
#-----------------------------------------#
#-----------------------------------------#

# Load run historic model file.
include(joinpath(@__DIR__, "..", "calibration", "run_historic_models", "run_dice_historic_climate_annual.jl"))

# Load log-posterior file for DICE climate model.
include(joinpath(@__DIR__, "..", "calibration", "create_log_posterior_dice_climate.jl"))

# Load inital parameter values for the DICE climate model (annualized version).
initial_parameters = DataFrame(load(joinpath(@__DIR__, "..", "data", "calibration_data", "calibration_initial_values_dice_climate.csv"), skiplines_begin=6))

# Create `DICE` function used in log-posterior calculations.
run_dice_climate! = construct_run_dice_climate_annual(calibration_end_year)

#----------------------------------#
#------ Baseline Calibration ------#
#----------------------------------#

!! Double check if any transition coefficients sum up to > 1 (would just be mid layer I think?) Could filter this out... or set this to -Inf in calibration

println("Begin calibration of annualized DICE climate model.\n")

# Create log-posterior function.
log_posterior_dice_climate = construct_dice_climate_log_posterior(run_dice_climate!, model_start_year=1850, calibration_end_year=calibration_end_year)

# Carry out Bayesian calibration using robust adaptive metropolis MCMC algorithm.
chain_dice_climate, accept_rate_dice_climate, cov_matrix_dice_climate = RAM_sample(log_posterior_dice_climate, initial_parameters.starting_point, Diagonal(initial_parameters.mcmc_step_size), n_mcmc_samples, opt_α=0.234)

# Discard burn-in values.
burned_chain_dice_climate = chain_dice_climate[Int(burn_in_length+1):end, :]

# Calculate mean posterior parameter values.
mean_dice_climate = vec(mean(burned_chain_dice_climate, dims=1))

# Create thinned chains (after burn-in period) with 10,000 and 100,000 samples and assign parameter names to each column.
thin100k_chain_dice_climate = DataFrame(burned_chain_dice_climate[thin_indices_100k, :])
thin10k_chain_dice_climate  = DataFrame(burned_chain_dice_climate[thin_indices_10k, :])

names!(thin100k_chain_dice_climate, [Symbol(initial_parameters.parameter[i]) for i in 1:length(mean_dice_climate)])
names!(thin10k_chain_dice_climate,  [Symbol(initial_parameters.parameter[i]) for i in 1:length(mean_dice_climate)])

!!! TO DO, SET UP OUTPUT FOLDER DIRECTORY FOR SAVING RESULTS !!!
output = joinpath("results", "dice_climate_params_may")
mkdir(joinpath(output, "calibrated_parameters"))
save(joinpath(@__DIR__, output, "calibrated_parameters", "mcmc_acceptance_rate.csv"), DataFrame(page_acceptance=accept_rate_dice_climate))
save(joinpath(@__DIR__, output, "calibrated_parameters", "mean_parameters.csv"), DataFrame(parameter = initial_parameters.parameter[1:length(mean_dice_climate)], mean=mean_dice_climate))
save(joinpath(@__DIR__, output, "calibrated_parameters", "parameters_10k.csv"), DataFrame(thin10k_chain_dice_climate))
save(joinpath(@__DIR__, output, "calibrated_parameters", "parameters_100k.csv"), DataFrame(thin100k_chain_dice_climate))



#--------------------------------------#
#--------------------------------------#
#------------ FUND Climate ------------#
#--------------------------------------#
#--------------------------------------#

# Load run historic model file.
include(joinpath(@__DIR__, "..", "calibration", "run_historic_models", "run_fund_historic_climate.jl"))

# Load log-posterior file for DICE climate model.
include(joinpath(@__DIR__, "..", "calibration", "create_log_posterior_fund_climate.jl"))

# Load inital parameter values for the DICE climate model (annualized version).
initial_parameters = DataFrame(load(joinpath(@__DIR__, "..", "data", "calibration_data", "calibration_initial_values_fund_climate.csv"), skiplines_begin=6))

# Create `DICE` function used in log-posterior calculations.
run_fund_climate! = construct_run_fund_climate(calibration_end_year)

#----------------------------------#
#------ Baseline Calibration ------#
#----------------------------------#

println("Begin calibration of FUND climate model.\n")

# Create log-posterior function.
log_posterior_fund_climate = construct_fund_climate_log_posterior(run_fund_climate!, model_start_year=1850, calibration_end_year=calibration_end_year)

# Carry out Bayesian calibration using robust adaptive metropolis MCMC algorithm.
chain_fund_climate, accept_rate_fund_climate, cov_matrix_fund_climate = RAM_sample(log_posterior_fund_climate, initial_parameters.starting_point, Diagonal(initial_parameters.mcmc_step_size), n_mcmc_samples, opt_α=0.234)

# Discard burn-in values.
burned_chain_fund_climate = chain_fund_climate[Int(burn_in_length+1):end, :]

# Calculate mean posterior parameter values.
mean_fund_climate = vec(mean(burned_chain_fund_climate, dims=1))

# Create thinned chains (after burn-in period) with 10,000 and 100,000 samples and assign parameter names to each column.
thin100k_chain_fund_climate = DataFrame(burned_chain_fund_climate[thin_indices_100k, :])
thin10k_chain_fund_climate  = DataFrame(burned_chain_fund_climate[thin_indices_10k, :])

names!(thin100k_chain_fund_climate, [Symbol(initial_parameters.parameter[i]) for i in 1:length(mean_fund_climate)])
names!(thin10k_chain_fund_climate,  [Symbol(initial_parameters.parameter[i]) for i in 1:length(mean_fund_climate)])

# !!! TO DO, SET UP OUTPUT FOLDER DIRECTORY FOR SAVING RESULTS !!!
output = joinpath("results", "dice_climate_params_may")
mkdir(joinpath(output, "calibrated_parameters"))
save(joinpath(@__DIR__, output, "calibrated_parameters", "mcmc_acceptance_rate.csv"), DataFrame(page_acceptance=accept_rate_fund_climate))
save(joinpath(@__DIR__, output, "calibrated_parameters", "mean_parameters.csv"), DataFrame(parameter = initial_parameters.parameter[1:length(mean_fund_climate)], mean=mean_fund_climate))
save(joinpath(@__DIR__, output, "calibrated_parameters", "parameters_10k.csv"), DataFrame(thin10k_chain_fund_climate))
save(joinpath(@__DIR__, output, "calibrated_parameters", "parameters_100k.csv"), DataFrame(thin100k_chain_fund_climate))



#----------------------------------------#
#----------------------------------------#
#------------ SNEASY + BRICK ------------#
#----------------------------------------#
#----------------------------------------#

# NOTE** This version uses the kernel density estimated marginal priors for the Antarctic ice sheet based on a calibration to paleo data.

# Load run historic model file.
include(joinpath(@__DIR__, "..", "calibration", "run_historic_models", "run_sneasy_brick_historic_climate.jl"))

# Load log-posterior file for SNEASY+BRICK model.
include(joinpath(@__DIR__, "..", "calibration", "create_log_posterior_sneasy_brick.jl"))

# Load inital parameter values for SNEASY+BRICK model.
initial_parameters = DataFrame(load(joinpath(@__DIR__, "..", "data", "calibration_data", "calibration_initial_values_sneasy_brick.csv"), skiplines_begin=6))

# Create `SNEASY+BRICK` function used in log-posterior calculations.
run_sneasybrick! = construct_run_sneasybrick(calibration_end_year)

#----------------------------------#
#------ Baseline Calibration ------#
#----------------------------------#

println("Begin calibration of SNEASY+BRICK model.\n")

# Create log-posterior function.
log_posterior_sneasy_brick = construct_sneasybrick_log_posterior(run_sneasybrick!, model_start_year=1850, calibration_end_year=calibration_end_year, joint_antarctic_prior=false)

# Carry out Bayesian calibration using robust adaptive metropolis MCMC algorithm.
chain_sneasybrick, accept_rate_sneasybrick, cov_matrix_sneasybrick = RAM_sample(log_posterior_sneasy_brick, initial_parameters.starting_point, Diagonal(initial_parameters.mcmc_step_size), n_mcmc_samples, opt_α=0.234)

# Discard burn-in values.
burned_chain_sneasybrick = chain_sneasybrick[Int(burn_in_length+1):end, :]

# Calculate mean posterior parameter values.
mean_sneasybrick = vec(mean(burned_chain_sneasybrick, dims=1))

# Create thinned chains (after burn-in period) with 10,000 and 100,000 samples and assign parameter names to each column.
thin100k_chain_sneasybrick = DataFrame(burned_chain_sneasybrick[thin_indices_100k, :])
thin10k_chain_sneasybrick  = DataFrame(burned_chain_sneasybrick[thin_indices_10k, :])

names!(thin100k_chain_sneasybrick, [Symbol(initial_parameters.parameter[i]) for i in 1:length(mean_sneasybrick)])
names!(thin10k_chain_sneasybrick,  [Symbol(initial_parameters.parameter[i]) for i in 1:length(mean_sneasybrick)])

# !!! TO DO, SET UP OUTPUT FOLDER DIRECTORY FOR SAVING RESULTS !!!
output = joinpath("results", "sneasybrick_params_may30")
mkdir(joinpath(output, "calibrated_parameters"))
save(joinpath(@__DIR__, output, "calibrated_parameters", "mcmc_acceptance_rate.csv"), DataFrame(page_acceptance=accept_rate_sneasybrick))
save(joinpath(@__DIR__, output, "calibrated_parameters", "mean_parameters.csv"), DataFrame(parameter = initial_parameters.parameter[1:length(mean_sneasybrick)], mean=mean_sneasybrick))
save(joinpath(@__DIR__, output, "calibrated_parameters", "parameters_10k.csv"), DataFrame(thin10k_chain_sneasybrick))
save(joinpath(@__DIR__, output, "calibrated_parameters", "parameters_100k.csv"), DataFrame(thin100k_chain_sneasybrick))

