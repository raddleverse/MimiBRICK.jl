# module TestCalibration

using MimiBRICK
using Test

import MimiBRICK: get_model, create_sneasy_brick, create_brick_doeclim
import MimiBRICK: construct_run_brick, construct_run_doeclimbrick, construct_run_sneasybrick
import MimiBRICK: construct_brick_log_posterior, construct_doeclimbrick_log_posterior, construct_sneasybrick_log_posterior
import MimiBRICK: run_calibration

##==============================================================================
## Checking short calibration, that it does things and isn't just perpetually stuck

calibration_start_year = 1850
calibration_end_year   = 2017
total_chain_length     = 1000
size_subsample         = 100
threshold_gr           = 1.1

# Create the log-posterior functions
log_posterior_brick = construct_brick_log_posterior(construct_run_brick(calibration_start_year, calibration_end_year), model_start_year=calibration_start_year, calibration_end_year=calibration_end_year, joint_antarctic_prior=false)
log_posterior_doeclimbrick = construct_doeclimbrick_log_posterior(construct_run_doeclimbrick(calibration_start_year, calibration_end_year), model_start_year=calibration_start_year, calibration_end_year=calibration_end_year, joint_antarctic_prior=false)
log_posterior_sneasybrick = construct_sneasybrick_log_posterior(construct_run_sneasybrick(calibration_start_year, calibration_end_year), model_start_year=calibration_start_year, calibration_end_year=calibration_end_year, joint_antarctic_prior=false)

# BRICK calibration
nparameters_brick = 35
x1 = run_calibration(log_posterior_brick; model_config="brick", calibration_start_year=1850, calibration_end_year=2017,
                     total_chain_length=total_chain_length, burnin_length=0, threshold_gr=threshold_gr, num_walkers=2,
                     size_subsample=size_subsample, start_from_priors=false)
@test size(x1[1])[1]==1000
@test size(x1[1])[2]==nparameters_brick
@test size(x1[5])[1]==100
@test size(x1[5])[2]==nparameters_brick
@test all(isa.(Matrix(x1[1]), Number))
@test !all(diff(x1[1][:,1]) .== 0)

# DOECLIM-BRICK calibration
nparameters_doeclimbrick = 44
x2 = run_calibration(log_posterior_doeclimbrick; model_config="doeclimbrick", calibration_start_year=1850, calibration_end_year=2017,
                     total_chain_length=total_chain_length, burnin_length=0, threshold_gr=threshold_gr, num_walkers=2,
                     size_subsample=size_subsample, start_from_priors=false)
@test size(x2[1])[1]==1000
@test size(x2[1])[2]==nparameters_doeclimbrick
@test size(x2[5])[1]==100
@test size(x2[5])[2]==nparameters_doeclimbrick
@test all(isa.(Matrix(x2[1]), Number))
@test !all(diff(x2[1][:,1]) .== 0)

# SNEASY-BRICK calibration
nparameters_sneasybrick = 51
x3 = run_calibration(log_posterior_sneasybrick; model_config="sneasybrick", calibration_start_year=1850, calibration_end_year=2017,
                     total_chain_length=total_chain_length, burnin_length=0, threshold_gr=threshold_gr, num_walkers=2,
                     size_subsample=size_subsample, start_from_priors=false)
@test size(x3[1])[1]==1000
@test size(x3[1])[2]==nparameters_sneasybrick
@test size(x3[5])[1]==100
@test size(x3[5])[2]==nparameters_sneasybrick
@test all(isa.(Matrix(x3[1]), Number))
@test !all(diff(x3[1][:,1]) .== 0)

end