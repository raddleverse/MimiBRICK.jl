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

output_dir = joinpath(@__DIR__, "results")
mkpath(output_dir)

calibration_start_year = 1850
calibration_end_year   = 2020
total_chain_length     = 15_000_000
burnin_length          = 5_000_000
num_walkers            = 2
size_subsample         = 10_000
threshold_gr           = 1.1
glacier_model          = :gsic # :gsic for Wigley-Raper; :mengel for Mengel et al
gmsl_data              = :wa # :wa for Wang et al; :cw for Church and White
glacier_data           = :fr # :dm for Dyurgerov-Meier; :fr for Frederikse et al. (2020)

# BRICK calibration
x = MimiBRICK.run_calibration(output_dir = output_dir, model_config="brick", calibration_start_year=calibration_start_year,
                    calibration_end_year=calibration_end_year, total_chain_length=total_chain_length,
                    burnin_length=burnin_length, threshold_gr=threshold_gr, num_walkers=num_walkers,
                    size_subsample=size_subsample, start_from_priors=false, joint_ais_prior=false, check_gr=true, 
                    glacier_model=:mengel, glacier_data=glacier_data, gmsl_data=gmsl_data)

# DOECLIM-BRICK calibration
x = MimiBRICK.run_calibration(output_dir = output_dir, model_config="doeclimbrick", calibration_start_year=calibration_start_year,
                    calibration_end_year=calibration_end_year, total_chain_length=total_chain_length,
                    burnin_length=burnin_length, threshold_gr=threshold_gr, num_walkers=num_walkers,
                    size_subsample=size_subsample, start_from_priors=false, joint_ais_prior=false, check_gr=true, 
                    glacier_model=:mengel, glacier_data=glacier_data, gmsl_data=gmsl_data)

# SNEASY-BRICK calibration
x = MimiBRICK.run_calibration(output_dir = output_dir, model_config="sneasybrick", calibration_start_year=calibration_start_year,
                    calibration_end_year=calibration_end_year, total_chain_length=total_chain_length,
                    burnin_length=burnin_length, threshold_gr=threshold_gr, num_walkers=num_walkers,
                    size_subsample=size_subsample, start_from_priors=false, joint_ais_prior=false, check_gr=true, 
                    glacier_model=:mengel, glacier_data=glacier_data, gmsl_data=gmsl_data)
