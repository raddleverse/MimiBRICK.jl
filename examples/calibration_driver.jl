using MimiBRICK

##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
## This file carries out a Markov chain Monte Carlo calibration of BRICK.
## This includes one of the following possible model configurations:
## (1) BRICK standalone (forced by input global mean surface temperatures and ocean heat uptake data)
## (2) DOECLIM+BRICK
## (3) SNEASY+BRICK
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------

output_dir = joinpath(@__DIR__, "..", "results")
mkpath(output_dir)

calibration_start_year = 1850
calibration_end_year   = 2017
total_chain_length     = 10_000_000 # 10M is sufficient for `brick` standalone; more for `doeclimbrick` or `sneasybrick`
burnin_length          = Int(0.2 * total_chain_length)
num_walkers            = 4
size_subsample         = 10_000 # good to check the ESS too (MCMCDiagnosticTools fcn `ess`)
threshold_gr           = 1.1

# BRICK calibration
x = MimiBRICK.run_calibration(output_dir = output_dir, model_config="brick", calibration_start_year=calibration_start_year,
                    calibration_end_year=calibration_end_year, total_chain_length=total_chain_length,
                    burnin_length=burnin_length, threshold_gr=threshold_gr, num_walkers=num_walkers,
                    size_subsample=size_subsample, start_from_priors=false)

# DOECLIM-BRICK calibration
x = MimiBRICK.run_calibration(output_dir = output_dir, model_config="doeclimbrick", calibration_start_year=calibration_start_year,
                    calibration_end_year=calibration_end_year, total_chain_length=total_chain_length,
                    burnin_length=burnin_length, threshold_gr=threshold_gr, num_walkers=num_walkers,
                    size_subsample=size_subsample, start_from_priors=false)

# SNEASY-BRICK calibration
x = MimiBRICK.run_calibration(output_dir = output_dir, model_config="sneasybrick", calibration_start_year=calibration_start_year,
                    calibration_end_year=calibration_end_year, total_chain_length=total_chain_length,
                    burnin_length=burnin_length, threshold_gr=threshold_gr, num_walkers=num_walkers,
                    size_subsample=size_subsample, start_from_priors=false)
