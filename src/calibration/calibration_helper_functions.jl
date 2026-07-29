using Missings
using DataFrames
using Distributions
using NetCDF
using KernelDensity
using CSVFiles
using Statistics
using XLSX
using CSV

#-------------------------------------------------------------------------------
# This file contains functions that are used for the model calibrations.
#-------------------------------------------------------------------------------

"""
    load_calibration_data(model_start_year::Int, last_calibration_year::Int; last_sea_level_norm_year::Int=1990, calibration_data_dir::Union{Nothing, String} = nothing, gmsl_data::Symbol=:cw)

Load and clean up data used for model calibration.

Description: This function loads, cleans up, and merges all of the calibration data into a single dataframe.

Function Arguments:

    - model_start_year         = The first year to include in the calibration data set.
    - last_calibration_year    = The last year to run the model for calibration (i.e. 1980 will not consider any post-1980 observations).
    - last_sea_level_norm_year = Some sea level data sets may need to be normalized to different years depending on when the calibration ends (this
                                 may be necessary for out-of-sample tests). These data sets will be normalized from 1961-last norm year, default = 1961-1990.
    - calibration_data_dir    = Data directory for calibration data. Defaults to package calibration data directory, changing this is not recommended.
    - gmsl_data               = GMSL reconstruction to use: `:wa` for Wang et al. (2024) or `:cw` for Church & White (2011).
    - calibration_targets     = `:standard` or `:mengel_ext`.
"""
function load_calibration_data(model_start_year::Int, last_calibration_year::Int; last_sea_level_norm_year::Int=1990, calibration_data_dir::Union{Nothing, String} = nothing, gmsl_data::Symbol=:cw, calibration_targets::Symbol=:standard)

    gmsl_data in (:wa, :cw) || throw(ArgumentError("gmsl_data must be :wa or :cw; got :$gmsl_data"))
    calibration_targets in (:standard, :mengel_ext) ||
        throw(ArgumentError("calibration_targets must be :standard or :mengel_ext; got :$calibration_targets"))

    # Create column of calibration years and calculate indicies for calibration time period relative to 1765-2020 (will crop later).
    # Note: first year is first year to run model (not necessarily year of first observation).
    df = DataFrame(year = collect(1765:max(2020, last_calibration_year)))
    model_calibration_indices = findall((in)(collect(model_start_year:last_calibration_year)), df.year)
    
    # set calibration data directory if one was not provided ie. it is set as nothing
    if isnothing(calibration_data_dir)
        calibration_data_dir = joinpath(@__DIR__, "..", "..", "data", "calibration_data")
    end    
    
    #-------------------------------------------------------------------
    # HadCRUT4 temperature data (anomalies relative to 1861-1880 mean).
    #-------------------------------------------------------------------

    # Load raw temperature data.
    raw_temp_data = DataFrame(load(joinpath(calibration_data_dir, "global_temp_hadcrut4.csv"), skiplines_begin=24))

    # Find indices to normalize temperature data to 1861-1880 mean.
    hadcrut_norm_indices = findall((in)(1861:1880), raw_temp_data[!,:year])

    # Normalize temperature data to 1861-1880 mean.
    norm_temp_data  = DataFrame(year=raw_temp_data[!,:year], hadcrut_temperature_obs = raw_temp_data[!,:median] .- mean(raw_temp_data[hadcrut_norm_indices, :median]))

    # Join data on year.
    #df = join(df, norm_temp_data, on=:year, kind=:outer)
    df = outerjoin(df, norm_temp_data, on=:year)

    # Read in HadCRUT4 1σ errors and rename column.
    raw_temp_errors  = DataFrame(load(joinpath(calibration_data_dir, "global_temp_hadcrut4_1sigma_uncertainty.csv"), skiplines_begin=21))
    rename!(raw_temp_errors, :one_sigma_all => :hadcrut_temperature_sigma)

    # Join data on year.
    #df = join(df, raw_temp_errors[!, [:year, :hadcrut_temperature_sigma]], on=:year, kind=:outer)
    df = outerjoin(df, raw_temp_errors[!, [:year, :hadcrut_temperature_sigma]], on=:year)

    #--------------------------------------------------------
    # Mauna Loa Instrumental Atmospheric CO₂ Concentrations.
    #--------------------------------------------------------

    # Load Mauna Loa CO₂ observations and errors, and rename columns.
    raw_mauna_loa_co2_data  = DataFrame(load(joinpath(calibration_data_dir, "co2_mauna_loa.csv"), skiplines_begin=59))
    rename!(raw_mauna_loa_co2_data, :mean => :maunaloa_co2_obs, :unc => :maunaloa_co2_sigma)

    # Join data on year.
    #df = join(df, raw_mauna_loa_co2_data, on=:year, kind=:outer)
    df = outerjoin(df, raw_mauna_loa_co2_data, on=:year)

    #-----------------------------------------------------
    # Law Dome Ice Core Atmospheric CO₂ Concentrations.
    #-----------------------------------------------------

    # Load Law Dome CO₂ observations and errors, and rename columns.
    raw_law_dome_co2_data = DataFrame(load(joinpath(calibration_data_dir, "law_dome_co2.csv"), skiplines_begin=4))
    rename!(raw_law_dome_co2_data, :co2_ice => :lawdome_co2_obs, :one_sigma_error => :lawdome_co2_sigma)

    # Join data on year.
    #df = join(df, raw_law_dome_co2_data, on=:year, kind=:outer)
    df = outerjoin(df, raw_law_dome_co2_data, on=:year)

    #---------------------------------------------------------------------------------
    # Annual Global Ocean Heat Content (0-3000 m).
    #---------------------------------------------------------------------------------

    # Load ocean heat content (0-3000m) observations and errors.
    ocean_heat_raw = DataFrame(load(joinpath(calibration_data_dir, "ocean_heat_gouretski_3000m.csv"), colnames=["year", "ocean_heat_obs", "ocean_heat_sigma"], skiplines_begin=3))

    # Join data on year.
    #df = join(df, ocean_heat_raw, on=:year, kind=:outer)
    df = outerjoin(df, ocean_heat_raw, on=:year)

    #---------------------------------------------------------------------------------
    # Decadal Ocean Carbon Fluxes
    #---------------------------------------------------------------------------------

    # Observations and errors from McNeil et al. (2003).
    ocean_co2_flux_data = DataFrame(year=[1985, 1995], oceanco2_flux_obs=[1.6, 2.0], oceanco2_flux_sigma=[0.4, 0.4])

    # Join data on year.
    #df = join(df, ocean_co2_flux_data, on=:year, kind=:outer)
    df = outerjoin(df, ocean_co2_flux_data, on=:year)

    #---------------------------------------------------------
    # Load the selected Global Mean Sea Level reconstruction.
    #---------------------------------------------------------

    if gmsl_data == :wa
        # Wang et al. (2024), an updated Church & White reconstruction spanning
        # 1900-2019. The source data are annual values and 1-sigma uncertainties
        # in millimeters (https://doi.org/10.5281/zenodo.15301866).
        wang_gmsl_file = joinpath(calibration_data_dir, "GMSL_yr.txt")
        if !isfile(wang_gmsl_file)
            wang_gmsl_url = "https://zenodo.org/records/15301866/files/GMSL_yr.txt?download=1"
            download(wang_gmsl_url, wang_gmsl_file)
        end

        # Use the original whitespace-delimited Zenodo file directly. Its three
        # columns are annual midpoint, GMSL (mm), and 1-sigma uncertainty (mm).
        raw_gmsl_data = DataFrame(CSV.File(
            wang_gmsl_file;
            header=[:year, :gmsl_mm, :gmsl_sigma_mm],
            delim=' ',
            ignorerepeated=true,
        ))
        gmsl_years = floor.(Int64, raw_gmsl_data.year)
        gmsl_obs   = raw_gmsl_data.gmsl_mm ./ 1000
        gmsl_error = raw_gmsl_data.gmsl_sigma_mm ./ 1000
    else
        # Church & White data updated through 2013.
        raw_gmsl_data = DataFrame(load(joinpath(calibration_data_dir, "CSIRO_Recons_gmsl_yr_2015.csv"), skiplines_begin=9))
        gmsl_obs   = raw_gmsl_data[!, Symbol("GMSL (mm)")] ./ 1000
        gmsl_error = raw_gmsl_data[!, Symbol("GMSL 1-sigma uncertainty (mm)")] ./ 1000
        # Dataset years are in half-years "so it is unambiguous as to which year they mean."
        gmsl_years = floor.(Int64, raw_gmsl_data[!, :Time])
    end

    # Get the indices for 1961-1990 and normalize GMSL observations relative to this period.
    gmsl_norm_indices = findall((in)(1961:last_sea_level_norm_year), gmsl_years)
    gmsl_obs_norm     = gmsl_obs .- mean(gmsl_obs[gmsl_norm_indices])

    # Combine into dataframe and join with other calibration data on year.
    gmsl_df = DataFrame(year = gmsl_years, gmsl_obs = gmsl_obs_norm, gmsl_sigma = gmsl_error)
    #df = join(df, gmsl_df, on=:year, kind=:outer)
    df = outerjoin(df, gmsl_df, on=:year)

    #---------------------------------------------------------------------------------
    # Load Greenland Ice Sheet Mass Balance Data for SIMPLE model.
    #---------------------------------------------------------------------------------

    # Load raw Greenland ice sheet mass balance data.
    raw_greenland_data = DataFrame(load(joinpath(calibration_data_dir, "greenland_MAR_InSAR_1958_2013.csv")))

    # Calculate number of observations (data set has different length montly and annual data in same sheet).
    n_obs = length(collect(skipmissing(raw_greenland_data[!, Symbol("Time (yearly)")])))

    # Years for SIMPLE data (covers 1958-2009 period).
    greenland_years = raw_greenland_data[1:n_obs, Symbol("Time (yearly)")]

    # Isolate observations of annual mass balance in meters of sea level equivalence.
    greenland_obs = raw_greenland_data[1:n_obs, Symbol("Yearly Cumulative MB in SLE (m)")]

    # Single observation error in meters (correpsonds to ± 30 Gt).
    greenland_error = repeat([raw_greenland_data[1, Symbol("Error in m")]], n_obs)

    # Get the indices of 1992-2001 (pooling with IMBIE data below, which spans 1992-2018, so normalize to 1992-2001 ten year period).
    greenland_norm_indices = findall((in)(1992:2001), greenland_years)

    # Normalize sea level observations to be relative to 1992-2001.
    greenland_obs_norm = greenland_obs .- mean(greenland_obs[greenland_norm_indices])

    # Combine into data frame and join with other calibration data.
  #  greenland_df = DataFrame(year = greenland_years, greenland_obs = greenland_obs_norm, greenland_sigma = greenland_error)
  #  df = join(df, greenland_df, on=:year, kind=:outer)

    #---------------------------------------------------------------------------------
    # Load IMBIE Greenland Ice Sheet Mass Balance Data (1992-2018).
    #---------------------------------------------------------------------------------

    # Load raw IMBIE Greenland ice sheet mass balance data.
    raw_imbie_greenland_data = DataFrame(load(joinpath(calibration_data_dir, "IMBIE_greenland_ice_sheet_1992_2018.csv"), skiplines_begin=5))

    # Isolate columns with required observations and assign shorter column names.
    raw_imbie_greenland_data = dropmissing(DataFrame(year = raw_imbie_greenland_data.Year,
                                                      greenland_obs=raw_imbie_greenland_data[!, Symbol("Cumulative ice sheet mass change (mm sea level)")],
                                                      greenland_sigma=raw_imbie_greenland_data[!, Symbol("Rate of ice sheet mass change uncertainty (mm sea level/yr)")]))

    # Get number of unique years with data (12 observations per year).
    years = floor.(raw_imbie_greenland_data.year)
    n_years = length(unique(years))

    # Assign vectors to hold annual values.
    annual_imbie_years = collect(years[1]:years[end])
    annual_imbie_obs   = zeros(n_years)
    annual_imbie_sigma = zeros(n_years)

    # Loop through each year, calculate annual value, and convert from mm to meters.
    for t = 1:n_years
        # Get indices for a unique year.
        year_indices = findall(x -> x == annual_imbie_years[t], years)
        # Calculate mean values for this year (corresponds to mid-year average).
        annual_imbie_obs[t]   = mean(raw_imbie_greenland_data.greenland_obs[year_indices]) / 1000
        # Convert annualized sigmas to corresponding monthly value and then add together for annual value (each monthly value has units of mm/yr)..
        annual_imbie_sigma[t] = sqrt(sum((raw_imbie_greenland_data.greenland_sigma[year_indices] ./ sqrt(12)).^2)) / 1000
    end

    # Normalize observations to 1992-2001 mean.
    greenland_imbie_norm_indices = findall((in)(1992:2001), annual_imbie_years)
    norm_annual_imbie_obs = annual_imbie_obs .- mean(annual_imbie_obs[greenland_imbie_norm_indices])

    # Merge into a dataframe and combine with other data.
    # Removing this to use the Frederikse et al 2020 data (below) instead
  #  annual_greenland_imbie_df = DataFrame(year=annual_imbie_years, greenland_imbie_obs=norm_annual_imbie_obs, greenland_imbie_sigma=annual_imbie_sigma)
  #  df = join(df, annual_greenland_imbie_df, on=:year, kind=:outer)

    #---------------------------------------------------------------------------------
    # Load Frederikse et al 2020 sea level contributor data (for Greenland in particular)
    #---------------------------------------------------------------------------------

    # Load raw Frederikse data for all GMSL contributors.
    raw_frederikse_data_xlsx = XLSX.readxlsx(joinpath(calibration_data_dir, "global_basin_timeseries.xlsx"))
    raw_frederikse_data = raw_frederikse_data_xlsx["Global"][:]
    column_names = raw_frederikse_data[1,:]
    column_names[1] = "Year"
    raw_frederikse_data = DataFrame(raw_frederikse_data[2:end,:], column_names)
    
    # Assign vectors to hold annual values.
    # /1000 to convert from mm to m
    annual_frederikse_years = raw_frederikse_data.Year
    annual_frederikse_obs   = raw_frederikse_data[!, Symbol("Greenland Ice Sheet [mean]")]/1000
    annual_frederikse_sigma = (raw_frederikse_data[!, "Greenland Ice Sheet [upper]"] - raw_frederikse_data[!, "Greenland Ice Sheet [lower]"])/4/1000

    # Normalize observations to 1992-2001 mean.
    greenland_frederikse_norm_indices = findall((in)(1961:last_sea_level_norm_year), annual_frederikse_years)
    norm_annual_frederikse_obs = annual_frederikse_obs .- mean(annual_frederikse_obs[greenland_frederikse_norm_indices])

    #-------------------------------------------------------------------------------------------
    # Combine Two Greenland Sets of Observations (for convenience) and normalize to common year.
    #-------------------------------------------------------------------------------------------

    # Get years that span first dataset (1958-2009) and IMBIE dataset (1992-2018).
    merged_years = collect(greenland_years[1]:annual_imbie_years[end])

    # Get indices for first dataset from 1958-1991, then use IMBIE dataset from 1992-2018.
    merged_indices = findall((in)(1958:1991), greenland_years)

    # Create merged observation values.
    merged_greenland_obs   = vcat(greenland_obs_norm[merged_indices], norm_annual_imbie_obs)
    merged_greenland_sigma = vcat(greenland_error[merged_indices], annual_imbie_sigma)

    # Normalize pooled observations to 1961-X mean (pooled data may be relative to the original Greenland data depending on cutoff year.. needed for out-of-sample tests).
    greenland_merged_norm_indices = findall((in)(1961:last_sea_level_norm_year), merged_years)
    norm_greenland_merged_obs     = merged_greenland_obs .- mean(merged_greenland_obs[greenland_merged_norm_indices])

    # Create dataframe of merged data and combine with other calibration data.
    # using the Frederikse et al data instead of IMBIE
    merged_greenland_df = DataFrame(year=annual_frederikse_years, merged_greenland_obs=norm_annual_frederikse_obs, merged_greenland_sigma=annual_frederikse_sigma)
    df = outerjoin(df, merged_greenland_df, on=:year)

    #---------------------------------------------------------------------------------
    # Load Glacier and Small Ice Caps (GSIC) Data
    #---------------------------------------------------------------------------------

    # Load GSIC raw data.
    raw_glaciers_data = DataFrame(load(joinpath(calibration_data_dir, "glacier_small_ice_caps_1961_2003.csv"), skiplines_begin=1))

    # Convert observations (cummulative GSIC melt contribution to sea level rise) from mm to meters.
    glaciers_obs = raw_glaciers_data[!, Symbol("contribution to sea level cumulative (mm)")] ./ 1000

    # Convert observation errors from mm to meters.
    glaciers_error = raw_glaciers_data[!, Symbol("standard dev. (mm/yr)")] ./ 1000

    # Years for GSIC data (covers 1961-2003 period).
    years = raw_glaciers_data[:, 1]

    # Get year indices and normalize data relative to 1961-1990 mean.
    glaciers_norm_indices = findall((in)(1961:last_sea_level_norm_year), years)
    glaciers_obs_norm = glaciers_obs .- mean(glaciers_obs[glaciers_norm_indices])

    # Combine into data frame and join with other calibration data.
    glaciers_df = DataFrame(year = years, glaciers_obs = glaciers_obs_norm, glaciers_sigma = glaciers_error)
    #df = join(df, glaciers_df, on=:year, kind=:outer)
    df = outerjoin(df, glaciers_df, on=:year)

    #---------------------------------------------------------------------------------
    # Load IMBIE Antarctic Ice Sheet Mass Balance Data (1992-2017).
    #---------------------------------------------------------------------------------

    # Load raw IMBIE Antarctic ice sheet mass balance data.
    raw_imbie_antarctic_data = DataFrame(load(joinpath(calibration_data_dir, "IMBIE_antarctic_ice_sheet_1992_2017.csv"), skiplines_begin=5))

    # Year 2017 does not have observations spanning full year, so drop 2017 partial values.
    drop_index = findlast(x -> x == 2016, floor.(raw_imbie_antarctic_data.Year))
    raw_imbie_antarctic_data = raw_imbie_antarctic_data[1:drop_index, :]

    # Give columns shorter names for convenience.
    #names!(raw_imbie_antarctic_data, Symbol.(["year", "cumm_obs_Gt", "cumm_sigma_Gt", "cumm_obs_mm", "cumm_sigma_mm"]))
    rename!(raw_imbie_antarctic_data, Symbol.(["year", "cumm_obs_Gt", "cumm_sigma_Gt", "cumm_obs_mm", "cumm_sigma_mm"]))

    # Dataset provides cummulative observation error (mm). Convert to annual observation error (mm/year) using same approach found in Greenland IMBIE dataset above.
    antarctic_sigma = zeros(size(raw_imbie_antarctic_data, 1))

    # Convert cummulative to annual observation errors for each monthly timestep.
    # Note* Like Greenland data above, values for each month expressed as mm/yr so need to account for 12 months.
    # Note** Last year has decreasing cummulative error (implying measurements that year had negative error? For these final periods, hold error constant at previous positive value).
    for t = 1:length(antarctic_sigma)
        if t == 1
            # Convert initial value.
            antarctic_sigma[t] = raw_imbie_antarctic_data.cumm_sigma_mm[t] * sqrt(12)
        else
            # Convert subsequent cummulative values.
            try
                antarctic_sigma[t] = sqrt(12*(raw_imbie_antarctic_data.cumm_sigma_mm[t]^2 - raw_imbie_antarctic_data.cumm_sigma_mm[t-1]^2))
            catch
                antarctic_sigma[t] = antarctic_sigma[t-1]
            end
        end
    end

    # Get number of unique years with data (12 observations per year).
    years   = floor.(raw_imbie_antarctic_data.year)
    n_years = length(unique(years))

    # Assign vectors to hold annual values.
    annual_imbie_years = collect(years[1]:years[end])
    annual_imbie_obs   = zeros(n_years)
    annual_imbie_sigma = zeros(n_years)

    # Loop through each year, calculate annual value, and convert from mm to meters.
    for t = 1:n_years
        # Get indices for a unique year.
        year_indices = findall(x -> x == annual_imbie_years[t], years)
        # Calculate mean values for this year (corresponds to mid-year average).
        annual_imbie_obs[t]   = mean(raw_imbie_antarctic_data.cumm_obs_mm[year_indices]) / 1000
        # Convert annualized sigmas to corresponding monthly value and then add together for annual value (each monthly value has units of mm/yr)..
        annual_imbie_sigma[t] = sqrt(sum((antarctic_sigma[year_indices] ./ sqrt(12)).^2)) / 1000
    end

    # Normalize observations to 1992-2001 mean.
    antarctic_imbie_norm_indices = findall((in)(1992:2001), annual_imbie_years)
    norm_annual_imbie_obs = annual_imbie_obs .- mean(annual_imbie_obs[antarctic_imbie_norm_indices])

    # Merge into a dataframe and combine with other data.
    annual_antarctic_imbie_df = DataFrame(year=annual_imbie_years, antarctic_imbie_obs=annual_imbie_obs, antarctic_imbie_sigma=annual_imbie_sigma)
    #df = join(df, annual_antarctic_imbie_df, on=:year, kind=:outer)
    df = outerjoin(df, annual_antarctic_imbie_df, on=:year)

    #---------------------------------------------------------------------------------
    # Load Original BRICK Antarctic Ice Sheet (AIS) Data.
    #---------------------------------------------------------------------------------
    # Comments from original BRICK code:
    #     Set up Shepherd et al 2012 window (a la Ruckert et al 2017)
    #     1992 to 2011 trend from Shepherd et al. 2012 is -71 +/- 53 Gt per yr
    #     We want the cumulative sea-level equivalent in meters for the year 2002
    #     Note the conversion: 360Gt = 1mm SLE

    # Values and equations taken from original R/fortran BRICK code.
    antarctic_inst_rate = abs(-71/360)/1000

    # "Using the midpoint of the 19-year interval"
    time_years = 2002-1992

    # Caluclate observed AIS value.
    obs_antarctic_inst = antarctic_inst_rate * time_years

    # Calculate AIS 1σ error.
    antarctic_inst_rate_error = abs(-53/360)/1000

    # Note on total error: sqrt(10) because 10 years of potentially accumulated error:
    #      Total error^2 = year₁ error² + year₂ error² + ... year₁₀ error² = 10*year X error².
    antarctic_inst_error = sqrt(time_years) * antarctic_inst_rate_error

    # Combine into data frame and join with other calibration data.
    # Note: Calibration sets year for this AIS observation as 2012.
    antarctic_inst_df = DataFrame(year = 2012, antarctic_inst_obs = obs_antarctic_inst, antarctic_inst_sigma = antarctic_inst_error)
    #df = join(df, antarctic_inst_df, on=:year, kind=:outer)
    df = outerjoin(df, antarctic_inst_df, on=:year)

    #---------------------------------------------------------------------------------
    # Marcus Sarofim extended calibration targets.
    #---------------------------------------------------------------------------------
    if calibration_targets == :mengel_ext
        targets_path = joinpath(calibration_data_dir, "calibration_targets_brick_mengel_ext.csv")
        isfile(targets_path) || error(
            "Marcus extended calibration targets not found: $targets_path"
        )

        targets = CSV.read(targets_path, DataFrame)
        required_columns = [
            :year,
            :ais, :ais_lo, :ais_hi,
            :gsic, :gsic_lo, :gsic_hi,
            :gis, :gis_lo, :gis_hi,
            :steric, :steric_lo, :steric_hi,
            :lws, :lws_lo, :lws_hi,
            :dang, :dang_sig, :dang_lo, :dang_hi,
        ]
        missing_columns = setdiff(required_columns, Symbol.(names(targets)))
        isempty(missing_columns) || error(
            "Marcus target file is missing required columns: $(join(missing_columns, ", "))"
        )

        # Marcus's table is in centimeters and is already referenced to the
        # 1995-2005 mean. Convert to meters to match MimiBRICK output. Reported
        # lower/upper bounds are 5-95% intervals; follow Marcus's conversion to
        # one-sigma errors and retain his 0.05 cm uncertainty floor.
        band_sigma(lo, hi) =
            ismissing(lo) || ismissing(hi) ? missing :
            max((hi - lo) / (2 * 1.645), 0.05) / 100
        cm_to_m(x) = ismissing(x) ? missing : x / 100

        mengel_ext_targets = DataFrame(
            year = Int.(targets.year),
            mengel_ext_ais_obs = cm_to_m.(targets.ais),
            mengel_ext_ais_sigma = band_sigma.(targets.ais_lo, targets.ais_hi),
            mengel_ext_glaciers_obs = cm_to_m.(targets.gsic),
            mengel_ext_glaciers_sigma = band_sigma.(targets.gsic_lo, targets.gsic_hi),
            mengel_ext_greenland_obs = cm_to_m.(targets.gis),
            mengel_ext_greenland_sigma = band_sigma.(targets.gis_lo, targets.gis_hi),
            mengel_ext_thermal_obs = cm_to_m.(targets.steric),
            mengel_ext_thermal_sigma = band_sigma.(targets.steric_lo, targets.steric_hi),
            mengel_ext_lws_obs = cm_to_m.(targets.lws),
            mengel_ext_lws_sigma = band_sigma.(targets.lws_lo, targets.lws_hi),
            mengel_ext_gmsl_obs = cm_to_m.(targets.dang),
            mengel_ext_gmsl_lower = cm_to_m.(targets.dang_lo),
            mengel_ext_gmsl_upper = cm_to_m.(targets.dang_hi),
            mengel_ext_gmsl_sigma = map(
                (dang_sigma, lws_lo, lws_hi) ->
                    ismissing(dang_sigma) || ismissing(lws_lo) || ismissing(lws_hi) ?
                    missing :
                    sqrt((dang_sigma / 100)^2 + band_sigma(lws_lo, lws_hi)^2),
                targets.dang_sig,
                targets.lws_lo,
                targets.lws_hi,
            ),
        )
        df = outerjoin(df, mengel_ext_targets, on=:year)
    end

    #---------------------------------------------------------------------------------
    # Finalize Joint Calibration Data Set.
    #---------------------------------------------------------------------------------

    # Sort all calibration data by year.
    sort!(df, :year)

    # Crop data to appropriate calibration years.
    df = df[model_calibration_indices, :]

    #---------------------------------------------------------------------------------
    # Load seperate trends dataset for Antarctic Ice Sheet.
    #---------------------------------------------------------------------------------
    model_years = collect(model_start_year:last_calibration_year)

    # From original BRICK code (DAIS module):
    #   Set up trends from IPCC AR5 for fitting AIS simulation.
    #   trends.ais = [Column 1=trend ; Column 2-3=90% window ; Column 4-5=beginning/ending year ;
    #                 Column 6-7=beginning/ending model indices for trend period]
    #   Note: IPCC trends are through 2010, but for the hindcast calibration, the
    #         forcing data go through 2009.
    #   Note: Each row is a different trend

    # Add a check for if the calibration end year occurs before the trends (for sensitivity tests).
    if last_calibration_year >= 2009

        ais_trends = Real[0.08  -0.1  0.27  1992  2001  findfirst(x -> x == 1992, model_years) findfirst(x -> x == 2001, model_years)
                          0.40  0.20  0.61  2002  2009  findfirst(x -> x == 2002, model_years) findfirst(x -> x == 2009, model_years)]

                        # Note: This 3rd trend was in original code, but not included here (check with Tony if it's double counting data.)
                        #0.27  0.16  0.38  1993  2009  findfirst(x -> x == 1993, model_years) findfirst(x -> x == 2009, model_years)]

        # Assign trends to dataframe and add column names.
        ais_trend_df = DataFrame(ais_trends, :auto)
        rename!(ais_trend_df, Symbol.(["Trend", "Lower_90_Percent", "Upper_90_Percent", "Start_Year", "End_Year", "Start_Model_Index", "End_Model_Index"]))

    else

        # Return -Inf if calibrating to a period not covering these trends.
        ais_trend_df = -Inf

    end

    #---------------------------------------------------------------------------------
    # Load seperate trends dataset for thermal expansion.
    #---------------------------------------------------------------------------------
    model_years = collect(model_start_year:last_calibration_year)

    # From original BRICK code (TE module):
    #   Set trends from IPCC AR5 Ch 13 (observational) to match:
    #   te_trends = [Column 1=trend ; Column 2-3=90% window ; Column 4-5=beginning/ending year ;
    #                Column 6-7=beginning/ending model indices for trend period]
    #   Note: IPCC trends are through 2010, but for the hindcast calibration, the forcing data go through 2009.

    # Add a check for if the calibration end year occurs before the trends (for sensitivity tests).
    if last_calibration_year >= 2009

        te_trends  = Real[0.8  0.5  1.1  1971  2009  findfirst(x -> x == 1971, model_years) findfirst(x -> x == 2009, model_years)
                          1.1  0.8  1.4  1993  2009  findfirst(x -> x == 1993, model_years) findfirst(x -> x == 2009, model_years)]

        # Assign trends to dataframe and add column names.
        te_trend_df = DataFrame(te_trends, :auto)
        #names!(te_trend_df, Symbol.(["Trend", "Lower_90_Percent", "Upper_90_Percent", "Start_Year", "End_Year", "Start_Model_Index", "End_Model_Index"]))
        rename!(te_trend_df, Symbol.(["Trend", "Lower_90_Percent", "Upper_90_Percent", "Start_Year", "End_Year", "Start_Model_Index", "End_Model_Index"]))

    else

        # Return -Inf if calibrating to a period not covering these trends.
        te_trend_df = -Inf

    end


    #--------------------------------------------------
    # Return all calibration observations and trends.
    #--------------------------------------------------
    return df, ais_trend_df, te_trend_df
end

"""
    mengel_extended_loglikelihood(...)

Calculate the five-series sea-level likelihood used for the `:mengel_ext`
calibration-target option. Each annual series uses its own non-missing years.
All modeled components are re-referenced to 1995-2005. The total-GMSL
comparison uses modeled AIS + glaciers + GIS + steric plus observed LWS, as in
Marcus Sarofim's `brick-fm` calibration.
"""
function mengel_extended_loglikelihood(
    calibration_data::DataFrame,
    model_years,
    modeled_glaciers,
    modeled_greenland,
    modeled_antarctic,
    modeled_thermal,
    σ_glaciers,
    ρ_glaciers,
    σ_greenland,
    ρ_greenland,
    σ_antarctic,
    ρ_antarctic,
    σ_thermal,
    ρ_thermal,
    σ_gmsl,
    ρ_gmsl,
)
    norm_indices = findall((in)(1995:2005), model_years)
    length(norm_indices) == 11 ||
        error("The :mengel_ext calibration requires model years covering 1995-2005.")

    rereference(x) = x .- mean(x[norm_indices])
    glaciers = rereference(modeled_glaciers)
    greenland = rereference(modeled_greenland)
    antarctic = rereference(modeled_antarctic)
    thermal = rereference(modeled_thermal)

    function series_loglikelihood(model, obs_column, sigma_column, σ_process, ρ)
        indices = findall(
            i -> !ismissing(calibration_data[i, obs_column]) &&
                 !ismissing(calibration_data[i, sigma_column]),
            1:nrow(calibration_data),
        )
        isempty(indices) && error("No usable observations found for $obs_column.")
        residual = Float64[
            calibration_data[i, obs_column] - model[i] for i in indices
        ]
        obs_sigma = Float64[calibration_data[i, sigma_column] for i in indices]
        return hetero_logl_ar1(residual, σ_process, ρ, obs_sigma)
    end

    llik = series_loglikelihood(
        glaciers, :mengel_ext_glaciers_obs, :mengel_ext_glaciers_sigma, σ_glaciers, ρ_glaciers
    )
    llik += series_loglikelihood(
        greenland, :mengel_ext_greenland_obs, :mengel_ext_greenland_sigma, σ_greenland, ρ_greenland
    )
    llik += series_loglikelihood(
        antarctic, :mengel_ext_ais_obs, :mengel_ext_ais_sigma, σ_antarctic, ρ_antarctic
    )
    llik += series_loglikelihood(
        thermal, :mengel_ext_thermal_obs, :mengel_ext_thermal_sigma, σ_thermal, ρ_thermal
    )

    gmsl_indices = findall(
        i -> !ismissing(calibration_data[i, :mengel_ext_gmsl_obs]) &&
             !ismissing(calibration_data[i, :mengel_ext_gmsl_sigma]) &&
             !ismissing(calibration_data[i, :mengel_ext_lws_obs]),
        1:nrow(calibration_data),
    )
    isempty(gmsl_indices) && error("No usable Marcus/Dangendorf GMSL observations found.")
    total = antarctic .+ glaciers .+ greenland .+ thermal
    gmsl_residual = Float64[
        calibration_data[i, :mengel_ext_gmsl_obs] -
        (total[i] + calibration_data[i, :mengel_ext_lws_obs])
        for i in gmsl_indices
    ]
    gmsl_sigma = Float64[calibration_data[i, :mengel_ext_gmsl_sigma] for i in gmsl_indices]
    llik += hetero_logl_ar1(gmsl_residual, σ_gmsl, ρ_gmsl, gmsl_sigma)

    # Additional Gaussian cumulative-change constraints used by Marcus.
    index_1992 = findfirst(==(1992), model_years)
    index_2017 = findfirst(==(2017), model_years)
    index_1961 = findfirst(==(1961), model_years)
    index_2003 = findfirst(==(2003), model_years)
    any(isnothing, (index_1992, index_2017, index_1961, index_2003)) &&
        error("The :mengel_ext calibration requires model years 1961, 1992, 2003, and 2017.")
    llik += logpdf(
        Normal(0.0072, 0.00156),
        antarctic[index_2017] - antarctic[index_1992],
    )
    llik += logpdf(
        Normal(0.02127, 0.00148),
        glaciers[index_2003] - glaciers[index_1961],
    )

    return llik
end

"""
    mengel_ext_thermal_process_logprior(σ_thermal, ρ_thermal)

Prior for the additional steric AR(1) discrepancy parameters used only by the
`:mengel_ext` target set.
"""
function mengel_ext_thermal_process_logprior(σ_thermal, ρ_thermal)
    return logpdf(Uniform(1e-10, 0.05), σ_thermal) +
           logpdf(truncated(Normal(0.8, 0.25), -1.0, 1.0), ρ_thermal)
end


"""
    hetero_logl_ar1(residuals, σ, ρ, ϵ)

Calculate AR(1) log-likelihood.

Description: This function calculates the AR(1) log-likelihood in terms of the data-model residuls when accounting for
             time-varying observation errors. It follows "The Effects of Time-Varying Observation Errors on Semi-Empirical
             Sea-Level Projections" (Ruckert et al., 2017) DOI 10.1007/s10584-016-1858-z.

Function Arguments:

      residuals = A vector of data-model residuals.
      σ         = AR(1) innovation standard deviation.
      ρ         = AR(1) autocorrelation term.
      ϵ         = A vector of time-varying observation error estimates (from calibration data sets).
"""
function hetero_logl_ar1(
    residuals::AbstractVector{<:Real},
    σ::Real,
    ρ::Real,
    ϵ::AbstractVector,
)

    # Calculate length of residuals.
    n=length(residuals)

    # Define AR(1) stationary process variance.
    σ_process = σ^2/(1-ρ^2)

    # Initialize AR(1) covariance matrix (just for convenience).
    H = abs.(collect(1:n)' .- collect(1:n))

    # Calculate residual covariance matrix (sum of AR(1) process variance and observation error variances).
    # Note: This follows Supplementary Information Equation (10) in Ruckert et al. (2017).
    cov_matrix = σ_process * ρ .^ H + Diagonal(Float64.(ϵ).^2)

    # Return the log-likelihood.
    return logpdf(MvNormal(cov_matrix), residuals)
end

"""

    hetero_logl_car1(residuals::Array{Float64,1}, indices::Array{Int64,1}, σ²_white_noise::Float64, α₀::Float64, ϵ::Array{Union{Float64, Missings.Missing},1})

Calculate CAR(1) log-likelihood.

Description: This function calculates the continuous time autoregressive, or CAR(1), log-likelihood for irregularly
             spaced data in terms of the data-model residuls when accounting for time-varying observation errors. It
             builds off of "The Effects of Time-Varying Observation Errors on Semi-Empirical Sea-Level Projections"
             (Ruckert et al., 2017) DOI 10.1007/s10584-016-1858-z and "The Analysis of Irregularly Observed Stochastic
             Astronomical Time-Series – I. Basics of Linear Stochastic Differential Equations" (Koen, 2005)
             doi.org/10.1111/j.1365-2966.2005.09213.x

Function Arguments:

      residuals      = A vector of data-model residuals.
      indices        = Index positions of observations relative to model time horizon (i.e. the first model time period = 1, the second = 2, etc.).
      σ²_white_noise = Variance of the continuous white noise process.
      α₀             = Parameter describing correlation memory of CAR(1) process.
      ϵ              = A vector of time-varying observation error estimates (from calibration data sets).
"""
function hetero_logl_car1(residuals::Array{Float64,1}, indices::Array{Int64,1}, σ²_white_noise::Float64, α₀::Float64, ϵ::Array{Union{Float64, Missings.Missing},1})

    # Calculate length of residuals.
    n=length(residuals)

    # Initialize covariance matrix for irregularly spaced data with relationships decaying exponentially.
    H = exp.(-α₀ .* abs.(indices' .- indices))

    # Define the variance of x(t), a continous stochastic time-series.
    σ² = σ²_white_noise / (2*α₀)

    # Calculate residual covariance matrix (sum of CAR(1) process variance and observation error variances).
    cov_matrix = σ² .* H + Diagonal(ϵ.^2)

    # Return the log-likelihood.
    return logpdf(MvNormal(cov_matrix), residuals)
end

"""
    calculate_trends(model_output::Array{Float64,1}, obs_trends::DataFrame, start_year::Int, end_year::Int)

Calculate modeled sea level trends.

Description: This function calculates the trend in different modeled sea level contributions over
             a user-specified time period.

Function Arguments:

      model_output = A vector of modeled annual sea level rise values.
      obs_trends   = The matrix of observed sea level trends provided by the "load_calibration_data" function defined above.
      start_year   = First year of the model run.
      end_year     = Last year of the model run.
"""
function calculate_trends(model_output::Array{Float64,1}, obs_trends::DataFrame, start_year::Int, end_year::Int)

    # Get number of trends to calculate.
    modeled_trends = zeros(size(obs_trends,1))

    # Loop through data for each different time period's trend.
    for i = 1:length(modeled_trends)

        # Find indices in model years corresponding to trend.
        trend_indices = findall((in)(obs_trends.Start_Year[i]:obs_trends.End_Year[i]), start_year:end_year)

        # "These calculate the least squares regression slope coefficients..." -original R/Fortran BRICK code.
        x = collect(trend_indices)
        y = model_output[trend_indices]

        # Calculate model trends (1000 × model trends because they're in meters, but observed trends in mm).
        modeled_trends[i] = 1000.0 * sum((x.-mean(x)) .* (y.-mean(y))) / sum((x.-mean(x)).^2)
    end
    return modeled_trends
end

"""
    truncated_kernel(data, lower_bound, upper_bound)

Calculate kernel density estimates with truncated bounds.

Description: This function creates a fitted kernel density object (from KernelDensity.jl package) and
             crops the edges to user-specified lower and upper bounds (setting boundaries beforehand,
             can lead to issues with some parameters because, "Due to the fourier transforms used
             internally, there should be sufficient spacing to prevent wrap-around at the boundaries")

Function Arguments:

      data        = Parameter samples used for the kernel density estimate.
      lower_bound = Minimum parameter value corresponding to the lower bound of the fitted density.
      upper_bound = Maximum parameter value corresponding to the upper bound of the fitted density.
"""
function truncated_kernel(data, lower_bound, upper_bound)

    # Estimate kernel density without specifying boundaries to avoid wrap-around at margins.
    kde_0 = kde(data)

    # Get number of values in KDE grid.
    n = length(kde_0.x)

    # Construct range of values that spans the KDE.
    kde_range = LinRange(minimum(kde_0.x), maximum(kde_0.x), n)

    # Get indices correponding to kde_range that represent desired truncated limits.
    lower_index = searchsortedfirst(kde_range, lower_bound)
    upper_index = searchsortedlast(kde_range, upper_bound)

    # Construct a truncated KDE object.
    truncated_kde = UnivariateKDE(kde_0.x[lower_index:upper_index], kde_0.density[lower_index:upper_index])

    # Interpolate the truncated KDE object for efficiency.
    # Note: "If you are making multiple calls to pdf, it will be more efficient to construct an intermediate InterpKDE to store the interpolation structure:"
    truncated_kde_interp = InterpKDE(truncated_kde)

    # Return truncated and interpolated KDE.
    return truncated_kde_interp
end

"""
    simulate_ar1_noise(n::Int, σ::Float64, ρ::Float64, ϵ::Array{Float64,1})

Simulate stationary AR(1) process with time varying observation errors.

Description: This function simulates a stationary AR(1) process (given time-varying observation errors supplied with
             each calibration data set) to superimpose noise onto the climate model projections.

Function Arguments:

      n = Number of time periods (years) the model is being run for.
      σ = Calibrated standard deviation.
      ρ = Calibrated autocorrelation coefficient.
      ϵ = Time-varying observation errors.
"""
function simulate_ar1_noise(n::Int, σ::Float64, ρ::Float64, ϵ::Array{Float64,1})

    # Define AR(1) stationary process variance.
    σ_process = σ^2/(1-ρ^2)

    # Initialize AR(1) covariance matrix (just for convenience).
    H = abs.(collect(1:n)' .- collect(1:n))

    # Calculate residual covariance matrix (sum of AR(1) process variance and observation error variances).
    # Note: This follows Supplementary Information Equation (10) in Ruckert et al. (2017).
    cov_matrix = σ_process * ρ .^ H + Diagonal(ϵ.^2)

    # Return a mean-zero AR(1) noise sample accounting for time-varying observation error.
    return rand(MvNormal(cov_matrix))
end

"""
    simulate_car1_noise(n, α₀, σ²_white_noise, ϵ)

Simulate stationary CAR(1) process with time varying observation errors.

Description: This function simulates a stationary CAR(1) process (given time-varying observation errors supplied with
             each calibration data set) to superimpose noise onto the climate model projections.

Function Arguments:

      n              = Number of time periods (years) the model is being run for.
      α₀             = Calibrated term describing correlation memory of CAR(1) process.
      σ²_white_noise = Calibrated continuous white noise process variance term.
      ϵ              = Time-varying observation errors.
"""
function simulate_car1_noise(n, α₀, σ²_white_noise, ϵ)

    # Indices for full time horizon.
    indices = collect(1:n)

    # Initialize covariance matrix for irregularly spaced data with relationships decaying exponentially.
    H = exp.(-α₀ .* abs.(indices' .- indices))

    # Define the variance of x(t), a continous stochastic time-series.
    σ² = σ²_white_noise / (2*α₀)

    # Calculate residual covariance matrix (sum of CAR(1) process variance and observation error variances).
    cov_matrix = σ² .* H + Diagonal(ϵ.^2)

    # Return a mean-zero CAR(1) noise sample accounting for time-varying observation error.
    return rand(MvNormal(cov_matrix))
end

"""
    replicate_errors(start_year::Int, end_year::Int, error_data)

REPLICATE TIME-VARYING OBSERVATION ERRORS FOR PERIODS WITHOUT COVERAGE.

Description: This function creates a time-series of observation errors for the entire model time horizon. For years
             without observation error estimates, the error remains constant at the average of the ten nearest error
             values in time.

Function Arguments:

      start_year = The first year to run the climate model.
      end_year   = The final year to run the climate model.
      error_data = A vector of time-varying observation errors supplied with each calibration data set.
"""
function replicate_errors(start_year::Int, end_year::Int, error_data)

    # Initialize model years and measurement error vector.
    model_years = collect(start_year:end_year)
    errors = zeros(length(model_years))

    # Find indices for periods that have observation errors.
    err_indices = findall(x-> !ismissing(x), error_data)

    # Replicate errors for all periods prior to start of observations based on average of 10 nearest errors in time.
    errors[1:(err_indices[1]-1)] .= mean(error_data[err_indices[1:10]])

    # Add errors for periods with observations.
    errors[err_indices[1:end]] = error_data[err_indices[1:end]]

    # Replicate errors for all periods after observation data based on average of 10 nearest errors in time.
    errors[(err_indices[end]+1):end] .= mean(error_data[err_indices[(end-9):end]])

    return errors
end

##------------------------------------------------------------------------------
## End
##------------------------------------------------------------------------------
