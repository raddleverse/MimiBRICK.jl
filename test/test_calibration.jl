@testitem "Calibration" begin
    # Checking short calibration, that it does things and isn't just perpetually stuck

    mktempdir() do tmp_dir
        calibration_start_year = 1850
        calibration_end_year   = 2020
        total_chain_length     = 1000
        size_subsample         = 100
        threshold_gr           = 1.1

        # BRICK calibration
        nparameters_brick = 35
        x1 = MimiBRICK.run_calibration(output_dir=tmp_dir, model_config="brick", calibration_start_year=calibration_start_year, calibration_end_year=calibration_end_year,
                            total_chain_length=total_chain_length, burnin_length=0, check_gr=false,
                            size_subsample=size_subsample, start_from_priors=false)
        @test size(x1[1])[1]==1000
        @test size(x1[1])[2]==nparameters_brick
        @test size(x1[5])[1]==100
        @test size(x1[5])[2]==nparameters_brick
        @test all(isa.(Matrix(x1[1]), Number))
        @test !all(diff(x1[1][:,1]) .== 0)
        @test isfile(joinpath(tmp_dir, "wigley-raper-glac", "parameters_subsample_brick_glac-ze_gmsl-wa.csv"))

        # DOECLIM-BRICK calibration
        nparameters_doeclimbrick = 44
        x2 = MimiBRICK.run_calibration(output_dir=tmp_dir, model_config="doeclimbrick", calibration_start_year=calibration_start_year, calibration_end_year=calibration_end_year,
                            total_chain_length=total_chain_length, burnin_length=0, check_gr=false,
                            size_subsample=size_subsample, start_from_priors=false)
        @test size(x2[1])[1]==1000
        @test size(x2[1])[2]==nparameters_doeclimbrick
        @test size(x2[5])[1]==100
        @test size(x2[5])[2]==nparameters_doeclimbrick
        @test all(isa.(Matrix(x2[1]), Number))
        @test !all(diff(x2[1][:,1]) .== 0)
        @test isfile(joinpath(tmp_dir, "wigley-raper-glac", "parameters_subsample_doeclimbrick_glac-ze_gmsl-wa.csv"))

        # SNEASY-BRICK calibration
        nparameters_sneasybrick = 51
        x3 = MimiBRICK.run_calibration(output_dir=tmp_dir, model_config="sneasybrick", calibration_start_year=calibration_start_year, calibration_end_year=calibration_end_year,
                            total_chain_length=total_chain_length, burnin_length=0, check_gr=false,
                            size_subsample=size_subsample, start_from_priors=false)
        @test size(x3[1])[1]==1000
        @test size(x3[1])[2]==nparameters_sneasybrick
        @test size(x3[5])[1]==100
        @test size(x3[5])[2]==nparameters_sneasybrick
        @test all(isa.(Matrix(x3[1]), Number))
        @test !all(diff(x3[1][:,1]) .== 0)
        @test isfile(joinpath(tmp_dir, "wigley-raper-glac", "parameters_subsample_sneasybrick_glac-ze_gmsl-wa.csv"))
    end
end

@testitem "Glacier calibration data options" begin
    using Statistics

    zemp_data, _, _ = MimiBRICK.load_calibration_data(1900, 2020)
    dm_data, _, _ = MimiBRICK.load_calibration_data(1900, 2020; glacier_data=:dm)

    zemp_indices = findall(x -> !ismissing(x), zemp_data.glaciers_obs)
    dm_indices = findall(x -> !ismissing(x), dm_data.glaciers_obs)

    @test zemp_data.year[zemp_indices] == collect(1962:2016)
    @test dm_data.year[dm_indices] == collect(1961:2003)
    @test length(unique(zemp_data.year[zemp_indices])) == length(zemp_indices)
    @test isapprox(mean(skipmissing(zemp_data[in.(zemp_data.year, Ref(1962:1990)), :glaciers_obs])), 0.0; atol=1e-12)
    @test all(>(0), zemp_data.glaciers_sigma[zemp_indices])
    @test zemp_data[zemp_data.year .== 2003, :glaciers_obs][1] != dm_data[dm_data.year .== 2003, :glaciers_obs][1]
    @test_throws ArgumentError MimiBRICK.load_calibration_data(1900, 2020; glacier_data=:invalid)
end

@testitem "GMSL calibration data options" begin
    using Statistics

    wang_data, _, _ = MimiBRICK.load_calibration_data(1900, 2019, gmsl_data=:wa)
    church_white_data, _, _ = MimiBRICK.load_calibration_data(1900, 2019, gmsl_data=:cw)
    wang_source_file = joinpath(dirname(pathof(MimiBRICK)), "..", "data", "calibration_data", "GMSL_yr.txt")

    wang_indices = findall(x -> !ismissing(x), wang_data.gmsl_obs)
    church_white_indices = findall(x -> !ismissing(x), church_white_data.gmsl_obs)

    @test wang_data.year[wang_indices] == collect(1900:2019)
    @test isfile(wang_source_file)
    @test church_white_data.year[church_white_indices[end]] == 2013
    norm_years = in.(wang_data.year, Ref(1961:1990))
    @test isapprox(mean(skipmissing(wang_data[norm_years, :gmsl_obs])), 0.0; atol=1e-12)
    @test isapprox(wang_data[wang_data.year .== 2019, :gmsl_sigma][1], 0.00965; atol=1e-12)
    @test wang_data[wang_data.year .== 1900, :gmsl_obs][1] != church_white_data[church_white_data.year .== 1900, :gmsl_obs][1]
    @test_throws ArgumentError MimiBRICK.load_calibration_data(1900, 2019, gmsl_data=:invalid)
end
