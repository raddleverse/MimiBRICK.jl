module TestDownscaling

# Test use cases for the downscaling function in MimiBRICK.jl.  The `src/downscale.jl` 
# file contains routines to downscale the BRICK global sea level projections to local. 
# This uses the sea-level "fingerprints" of [Slangen et al. (2014)](https://link.springer.com/article/10.1007/s10584-014-1080-9). 

using MimiBRICK
using Test

# set up directory 
tmp_dir=joinpath(@__DIR__, "tmp_dir")
mkpath(tmp_dir)

# calibration set up - very small set since we just need files to run downscaling
calibration_start_year  = 1850
calibration_end_year    = 2017
total_chain_length      = 100
size_subsample          = 10
threshold_gr            = 1.1

for model_config in ["brick", "sneasybrick", "doeclimbrick"]
    
    # run calibration
    MimiBRICK.run_calibration(output_dir=tmp_dir, model_config=model_config, calibration_start_year=1850, calibration_end_year=2017,
                        total_chain_length=total_chain_length, burnin_length=0, threshold_gr=threshold_gr, num_walkers=2,
                        size_subsample=size_subsample, start_from_priors=false)

    # run hindcast and projections 
    MimiBRICK.run_hindcast(output_dir=tmp_dir, model_config=model_config)
    MimiBRICK.run_projections(output_dir=tmp_dir, model_config=model_config)

    # Lat and Lon for New York City
    lat=40.7128 # deg N
    lon=360-74.0060 # 74.0060 deg W

    # specification 1
    # testing hindcast ensemble
    years, lsl_hind_ens=MimiBRICK.downscale_brick(lon=lon, lat=lat, results_dir=tmp_dir, proj_or_hind="hind", ensemble_or_map="ensemble", model_config=model_config, rcp_scenario="RCP26")
    @test size(lsl_hind_ens)[1]==168
    @test size(lsl_hind_ens)[2]==10
    @test length(years)==168
    @test all([isa(lsl_hind_ens[i,end],Number) for i=1:length(years)])
    @test all([isa(years[i],Number) for i=1:length(years)])

    # testing projections ensemble
    years, lsl_proj_ens=MimiBRICK.downscale_brick(lon=lon, lat=lat, results_dir=tmp_dir, proj_or_hind="proj", ensemble_or_map="ensemble", model_config=model_config, rcp_scenario="RCP85")
    @test size(lsl_proj_ens)[1]==451
    @test size(lsl_proj_ens)[2]==10
    @test length(years)==451
    @test all([isa(lsl_proj_ens[i,end],Number) for i=1:length(years)])
    @test all([isa(years[i],Number) for i=1:length(years)])

    # testing projections with MAP
    years, lsl_map=MimiBRICK.downscale_brick(lon=lon, lat=lat, results_dir=tmp_dir, proj_or_hind="proj", ensemble_or_map="map", model_config=model_config, rcp_scenario="RCP85")
    @test ndims(lsl_map)==1
    @test length(lsl_map)==451
    @test length(years)==451
    @test all([isa(lsl_map[i],Number) for i=1:length(years)])
    @test all([isa(years[i],Number) for i=1:length(years)])

end

# delete directory
rm(tmp_dir, recursive=true)

end
