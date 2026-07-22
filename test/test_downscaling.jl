@testitem "Downscaling" begin
    # Test use cases for the downscaling function in MimiBRICK.jl.  The `src/downscale.jl` 
    # file contains routines to downscale the BRICK global sea level projections to local. 
    # This uses the sea-level "fingerprints" of [Slangen et al. (2014)](https://link.springer.com/article/10.1007/s10584-014-1080-9). 


    # set up directory 
    mktempdir() do tmp_dir
        # calibration set up - very small set since we just need files to run downscaling
        calibration_start_year  = 1850
        calibration_end_year    = 2017
        total_chain_length      = 250
        size_subsample          = 10
        glacier_model           = :mengel

        for model_config in ["brick", "sneasybrick", "doeclimbrick"]
            
            # run calibration
            MimiBRICK.run_calibration(output_dir=tmp_dir, model_config=model_config, calibration_start_year=1850, calibration_end_year=2017,
                                total_chain_length=total_chain_length, burnin_length=0, check_gr=false,
                                size_subsample=size_subsample, start_from_priors=false, glacier_model=glacier_model)

            # run hindcast and projections 
            MimiBRICK.run_projections(output_dir=tmp_dir, model_config=model_config, ssprcp_scenario="ssp245", glacier_model=glacier_model)

            # Lat and Lon for New York City
            lat=40.7128 # deg N
            lon=360-74.0060 # 74.0060 deg W

            # testing projections ensemble
            years, lsl_proj_ens=MimiBRICK.downscale_brick(lon=lon, lat=lat, results_dir=tmp_dir, ensemble_or_map="ensemble", 
                                                          model_config=model_config, ssprcp_scenario="ssp245", glacier_model=glacier_model)
            @test size(lsl_proj_ens)[1]==451
            @test size(lsl_proj_ens)[2]==10
            @test length(years)==451
            @test all([isa(lsl_proj_ens[i,end],Number) for i=1:length(years)])
            @test all([isa(years[i],Number) for i=1:length(years)])

            # testing projections with MAP
            years, lsl_map=MimiBRICK.downscale_brick(lon=lon, lat=lat, results_dir=tmp_dir, ensemble_or_map="map", 
                                                     model_config=model_config, ssprcp_scenario="ssp245", glacier_model=glacier_model)
            @test ndims(lsl_map)==1
            @test length(lsl_map)==451
            @test length(years)==451
            @test all([isa(lsl_map[i],Number) for i=1:length(years)])
            @test all([isa(years[i],Number) for i=1:length(years)])

        end
    end
end
