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
#total_chain_length     = 20_000_000
#size_subsample         = 10_000
total_chain_length     = 200000
size_subsample         = 10
threshold_gr           = 1.1

##TESTING##
model_config="brick"; burnin_length=100; num_walkers=2; start_from_priors=false; calibration_data_dir = nothing;
#
x = MimiBRICK.run_calibration(output_dir = output_dir, model_config="brick", calibration_start_year=calibration_start_year,
                    calibration_end_year=calibration_end_year, total_chain_length=total_chain_length,
                    burnin_length=burnin_length, threshold_gr=threshold_gr, num_walkers=2,
                    size_subsample=size_subsample, start_from_priors=false)
##END TESTING##

# BRICK calibration
x = MimiBRICK.run_calibration(output_dir = output_dir, model_config="brick", calibration_start_year=calibration_start_year,
                    calibration_end_year=calibration_end_year, total_chain_length=total_chain_length,
                    burnin_length=1_000_000, threshold_gr=threshold_gr, num_walkers=2,
                    size_subsample=size_subsample, start_from_priors=false)

# DOECLIM-BRICK calibration
x = MimiBRICK.run_calibration(output_dir = output_dir, model_config="doeclimbrick", calibration_start_year=1850, calibration_end_year=2017,
                    total_chain_length=20_000_000, burnin_length=7_000_000, threshold_gr=1.1, num_walkers=2,
                    size_subsample=10_000, start_from_priors=false)

# SNEASY-BRICK calibration
x = MimiBRICK.run_calibration(output_dir = output_dir, model_config="sneasybrick", calibration_start_year=1850, calibration_end_year=2017,
                    total_chain_length=20_000_000, burnin_length=1_000_000, threshold_gr=1.1, num_walkers=2,
                    size_subsample=10_000, start_from_priors=false)
