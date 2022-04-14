
using Query
using CSV
using CSVFiles
using DataFrames
using NetCDF
using StatsBase

##==============================================================================
## Supporting Functions to Downscale BRICK from GMSL to LSL

"""
    get_fingerprints()
Retrieve BRICK fingerprints from NetCDF file - will download the file to a
folder `data` directory
"""
function get_fingerprints()

    fp_dir = joinpath(@__DIR__, "..", "data", "model_data")
    isdir(fp_dir) || mkpath(fp_dir)
    fp_file = joinpath(fp_dir, "FINGERPRINTS_SLANGEN_Bakker.nc")
    if !isfile(fp_file)
        url = "https://github.com/scrim-network/BRICK/raw/master/fingerprints/FINGERPRINTS_SLANGEN_Bakker.nc"
        download(url, fp_file)
    end

    fplat = ncread(fp_file,"lat")
    fplon = ncread(fp_file,"lon")
    fpAIS = ncread(fp_file,"AIS")
    fpGSIC = ncread(fp_file,"GLAC")
    fpGIS = ncread(fp_file,"GIS")
    ncclose()

    return fplat,fplon,fpAIS,fpGSIC,fpGIS
end

##==============================================================================
## Small Helper Functions for dealing with sea level fingerprints near land

"""
    next_lat(lat::Float64, inc::Int64, direction::Symbol)
Increment latitude by `inc` in either positive direction (`direction=:increase`)
or in the negative direction (`direction=:decrease`).
Assumes latitude runs from -90 to 90 (deg N).
"""
function next_lat(lat::Float64, inc::Int64, direction::Symbol)
    if lat < -90 || lat > 90
        error("Latitude must be between -90 and 90")
    end

    if direction == :increase
        new_lat = lat + inc
        if new_lat > 90
            new_lat = new_lat - 180 #wrap around
        end

    elseif direction == :decrease
        new_lat = lat - inc
        if new_lat < -90
            new_lat = new_lat + 180
        end
    end
    return new_lat
end

"""
    next_lon(lon::Float64, inc::Int64, direction::Symbol)
Increment longitude by `inc` in either positive direction
(`direction=:increase`) or in the negative direction (`direction=:decrease`).
Assumes longitude runs from 0 to 360 (deg E).
"""
function next_lon(lon::Float64, inc::Int64, direction::Symbol)
    if lon < 0 || lon > 360
        error("Longitude must be between 0 and 360")
    end

    if direction == :increase
        new_lon = lon + inc
        if new_lon > 360
            new_lon = new_lon - 360
        end
    elseif direction == :decrease
        new_lon = lon - inc
        if new_lon < 0
            new_lon = new_lon + 360
        end
    end
    return new_lon
end

##==============================================================================
"""
    downscale_brick(;lon, lat, proj_or_hind, ensemble_or_map, model_config, rcp_scenario="RCP85")
Downscale BRICK projections to a single point, using either the whole ensemble
or only the maximum a posteriori ensemble member
lon = longitude (degrees East) of location for downscaling
lat = latitude (degrees North) of location for downscaling
proj_or_hind = "proj" for projections, or "hind" for hindcast
ensemble_or_map = "ensemble" for entire posterior ensemble, or "map" for the maximum a posteriori ensemble member (single simulation)
model_config = "brick", "doeclimbrick", or "sneasybrick"
rcp_scenario = "RCP26", "RCP45", "RCP60", or "RCP85" (default). Doesn't matter for hindcast.
"""
function downscale_brick(;lon, lat, proj_or_hind, ensemble_or_map, model_config, rcp_scenario="RCP85")

    brick_results_dir = "my_brick_results_20M_20-02-2022"
    doeclimbrick_results_dir = "my_doeclimbrick_results_20M_19-02-2022"
    sneasybrick_results_dir = "my_sneasybrick_results_20M_19-02-2022"

    if model_config=="brick"
        results_dir = brick_results_dir
    elseif model_config=="doeclimbrick"
        results_dir = doeclimbrick_results_dir
    elseif model_config=="sneasybrick"
        results_dir = sneasybrick_results_dir
    end
    appen = "$(model_config)_$(results_dir[(findfirst("results_", results_dir)[end]+1):length(results_dir)])"

    if proj_or_hind=="proj"
        results_dir = joinpath(results_dir, "projections_csv", rcp_scenario)
        slr_dir = joinpath(@__DIR__, "..", "results", results_dir)
        MAP = DataFrame(load(joinpath(slr_dir,"projections_MAP_$(rcp_scenario)_$(appen).csv")))
        years = MAP[:,:YEAR]
        if ensemble_or_map=="ensemble"
            AIS = CSV.read(joinpath(slr_dir,"projections_antarctic_$(rcp_scenario)_$(appen).csv"), DataFrame)
            GIS = CSV.read(joinpath(slr_dir,"projections_greenland_$(rcp_scenario)_$(appen).csv"), DataFrame)
            GSIC = CSV.read(joinpath(slr_dir,"projections_glaciers_$(rcp_scenario)_$(appen).csv"), DataFrame)
            TE = CSV.read(joinpath(slr_dir,"projections_thermal_$(rcp_scenario)_$(appen).csv"), DataFrame)
            LWS = CSV.read(joinpath(slr_dir,"projections_landwater_storage_sl_$(rcp_scenario)_$(c).csv"), DataFrame)
            num_ens = size(AIS)[2]
        elseif ensemble_or_map=="map"
            AIS = MAP[:,:AIS]
            GIS = MAP[:,:GIS]
            GSIC = MAP[:,:GLAC]
            TE = MAP[:,:TE]
            LWS = MAP[:,:LWS]
            num_ens = 1
        end
    elseif proj_or_hind=="hind"
        results_dir = joinpath(results_dir, "hindcast_csv")
        slr_dir = joinpath(@__DIR__, "..", "results", results_dir)
        MAP = DataFrame(load(joinpath(slr_dir,"hindcast_MAP_$(appen).csv")))
        years = MAP[:,:YEAR]
        if ensemble_or_map=="ensemble"
            AIS = CSV.read(joinpath(slr_dir,"hindcast_antarctic_$(appen).csv"), DataFrame)
            GIS = CSV.read(joinpath(slr_dir,"hindcast_greenland_$(appen).csv"), DataFrame)
            GSIC = CSV.read(joinpath(slr_dir,"hindcast_glaciers_$(appen).csv"), DataFrame)
            TE = CSV.read(joinpath(slr_dir,"hindcast_thermal_$(appen).csv"), DataFrame)
            LWS = CSV.read(joinpath(slr_dir,"hindcast_landwater_storage_sl_$(appen).csv"), DataFrame)
            num_ens = size(AIS)[2]
        elseif ensemble_or_map=="map"
            AIS = MAP[:,:AIS]
            GIS = MAP[:,:GIS]
            GSIC = MAP[:,:GLAC]
            TE = MAP[:,:TE]
            LWS = MAP[:,:LWS]
            num_ens = 1
        end
    end
    num_years = length(years)

    (fplat,fplon,fpAIS,fpGSIC,fpGIS) = get_fingerprints()

    # Convert Longitude to degrees East
    # CIAM Lat is already in (-90,90) by default
    if lon <0
        lon = lon + 360
    end

    # Find fingerprint degrees nearest to lat,lon
    ilat = findall(isequal(minimum(abs.(fplat.-lat))),abs.(fplat.-lat))
    ilon = findall(isequal(minimum(abs.(fplon.-lon))),abs.(fplon.-lon))

    # Take average of closest lat/lon values
    fpAIS_flat = collect(skipmissing(Iterators.flatten(fpAIS[ilon,ilat])))
    fpGSIC_flat = collect(skipmissing(Iterators.flatten(fpGSIC[ilon,ilat])))
    fpGIS_flat = collect(skipmissing(Iterators.flatten(fpGIS[ilon,ilat])))

    fpAIS_loc = mean(fpAIS_flat[isnan.(fpAIS_flat).==false],dims=1)[1]
    fpGSIC_loc = mean(fpGSIC_flat[isnan.(fpGSIC_flat).==false],dims=1)[1]
    fpGIS_loc = mean(fpGIS_flat[isnan.(fpGIS_flat).==false],dims=1)[1]
    fpTE_loc = 1.0
    fpLWS_loc = 1.0

    # Keep searching nearby lat/lon values if fingerprint value is NaN unless limit is hit
    inc = 1

    while isnan(fpAIS_loc) || isnan(fpGIS_loc) || isnan(fpGSIC_loc) && inc<5

        newlonStart = next_lon.(fplon[ilon], inc, :decrease)[1]
        newlatStart = next_lat.(fplat[ilat], inc, :decrease)[1]
        newlonEnd = next_lon.(fplon[ilon], inc, :increase)[1]
        newlatEnd = next_lat.(fplat[ilat], inc, :increase)[1]

        latInd1 = minimum(findall(isequal(minimum(abs.(fplat.-newlatStart))),abs.(fplat.-newlatStart)))
        latInd2 = maximum(findall(isequal(minimum(abs.(fplat.-newlatEnd))),abs.(fplat.-newlatEnd)))

        lonInd1 = minimum(findall(isequal(minimum(abs.(fplon.-newlonStart))),abs.(fplon.-newlonStart)))
        lonInd2 = maximum(findall(isequal(minimum(abs.(fplon.-newlonEnd))),abs.(fplon.-newlonEnd)))

        if latInd2 < latInd1
            latInds = [latInd1; 1:latInd2]
        else
            latInds = latInd1:latInd2
        end

        if lonInd2 < lonInd1
            lonInds = [lonInd1; 1:lonInd2]
        else
            lonInds = lonInd1:lonInd2
        end

        fpAIS_flat = collect(skipmissing(Iterators.flatten(fpAIS[lonInds,latInds])))
        fpGSIC_flat = collect(skipmissing(Iterators.flatten(fpGSIC[lonInds,latInds])))
        fpGIS_flat = collect(skipmissing(Iterators.flatten(fpGIS[lonInds,latInds])))

        fpAIS_loc = mean(fpAIS_flat[isnan.(fpAIS_flat).==false],dims=1)[1]
        fpGSIC_loc = mean(fpGSIC_flat[isnan.(fpGSIC_flat).==false],dims=1)[1]
        fpGIS_loc = mean(fpGIS_flat[isnan.(fpGIS_flat).==false],dims=1)[1]

        inc = inc + 1

    end

    # If still NaN, throw an error
    if isnan(fpAIS_loc) || isnan(fpGIS_loc) || isnan(fpGSIC_loc)
        println("Error: no fingerprints found for (lon,lat) = ($(lon),$(lat))")
        return nothing
    end

    # Prepare outputs for CSV
    outputs = zeros(Union{Missing, Float64}, num_years, num_ens+1)
    outputs[:,1] = years

    # Multiply fingerprints by BRICK ensemble members
    if ndims(AIS) > 1
        lsl_out = zeros(Union{Missing, Float64}, size(AIS)[1], size(AIS)[2])
        for n in 1:size(AIS)[2] # loop through ensemble members
            lsl_out[:,n] = fpGIS_loc * GIS[:,n] + fpAIS_loc * AIS[:,n] + fpGSIC_loc * GSIC[:,n] +
                           fpTE_loc * TE[:,n] + fpLWS_loc * LWS[:,n]
        end
        outputs[:,2:end] = lsl_out
    else
        lsl_out = zeros(Union{Missing, Float64}, size(AIS)[1])
        lsl_out[:] = fpGIS_loc * GIS[:] + fpAIS_loc * AIS[:] + fpGSIC_loc * GSIC[:] +
                     fpTE_loc * TE[:] + fpLWS_loc * LWS[:]
        outputs[:,2] = lsl_out
    end

    # Write to CSV
    # NB: `my_brick_results_20M_20-02-2022` would need to be swapped out with whatever results directory your new runs are in, if you use a different ensemble
    filepath_output = joinpath(slr_dir, "localslr")
    mkpath(filepath_output)

    lat_rounded = round(lat, digits=2)
    lon_rounded = round(lon, digits=2)
    filename_output = joinpath(filepath_output,"projections_lsl-lat$(lat_rounded)-lon$(lon_rounded)_$(model_config).csv")
    CSV.write(filename_output, DataFrame(outputs, :auto))

    return years, lsl_out
end
