module TestFullRun

# These tests run full calibration for all three model configurations including
# calibration, hindcast and projections runs, and downscaling.  There are no explicit
# tests and instead just a check for clean runs without errors.

using MimiBRICK
using Test

# set up directory 
tmp_dir=joinpath(@__DIR__, "tmp_dir")
mkpath(tmp_dir)

# calibration set up
calibration_start_year  = 1850
calibration_end_year    = 2017
total_chain_length      = 1000
size_subsample          = 100
threshold_gr            = 1.1

for model_config in ["brick", "sneasybrick", "doeclimbrick"]

    # Create the log-posterior functions
    log_posterior_brick=MimiBRICK.construct_brick_log_posterior(MimiBRICK.construct_run_brick(calibration_start_year, calibration_end_year), model_start_year=calibration_start_year, calibration_end_year=calibration_end_year, joint_antarctic_prior=false)

    # run calibration
    MimiBRICK.run_calibration(log_posterior_brick; output_dir=tmp_dir, model_config=model_config, calibration_start_year=1850, calibration_end_year=2017,
                        total_chain_length=total_chain_length, burnin_length=0, threshold_gr=threshold_gr, num_walkers=2,
                        size_subsample=size_subsample, start_from_priors=false)

    # run hindcast and projections 
    MimiBRICK.run_hindcast(output_dir=tmp_dir, model_config=model_config)
    MimiBRICK.run_projections(output_dir=tmp_dir, model_config=model_config)

    # downscale

    # Lat and Lon for New York City
    lat=40.7128 # deg N
    lon=360-74.0060 # 74.0060 deg W

    years, lsl_hind_ens=MimiBRICK.downscale_brick(lon=lon, lat=lat, results_dir=tmp_dir, proj_or_hind="hind", ensemble_or_map="ensemble", model_config="brick", rcp_scenario="RCP85")
    years, lsl_proj_ens=MimiBRICK.downscale_brick(lon=lon, lat=lat, results_dir=tmp_dir, proj_or_hind="proj", ensemble_or_map="ensemble", model_config="brick", rcp_scenario="RCP85")
    years, lsl_map=MimiBRICK.downscale_brick(lon=lon, lat=lat, results_dir=tmp_dir, proj_or_hind="proj", ensemble_or_map="map", model_config="brick", rcp_scenario="RCP85")

end

# delete directory
rm(tmp_dir, recursive=true)

end # module
