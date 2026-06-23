using MimiBRICK
using CSV
using DataFrames
using MimiSNEASY
using Dates

##==============================================================================
## Read output from SNEASY RCP projections

dir_sneasybrick = joinpath(@__DIR__, "results", "projections_csv", "sneasybrick")

for scenario in ["ssp119","ssp126","ssp245","ssp370","ssp460","ssp585","ssp534-over"]

    filename_sneasy_map = joinpath(dir_sneasybrick,scenario,"projections_MAP_$(scenario)_sneasybrick.csv")
    map_projections = DataFrame(load(filename_sneasy_map))

    ## Get temperature and ocean heat
    years = map_projections[!,:YEAR]
    temperature = map_projections[!,:TEMP]
    ocean_heat = map_projections[!,:OCHEAT]

    ## Normalize temperature relative to 1850-1870
    ibeg = findall(x->x==1850,years)[1]
    iend = findall(x->x==1870,years)[1]
    temperature_norm = temperature .- mean(temperature[ibeg:iend])

    ## Write CSV files
    # temperature
    filename_temp = joinpath(@__DIR__, "data", "model_data", "sneasy_temperature_$(scenario)_$(Int(years[1]))_$(Int(years[end]))_$(Dates.format(now(),"dd-mm-yyyy")).csv")
    CSV.write(filename_temp, DataFrame([years,temperature_norm], ["Year","MAP Temperature"]))
    # ocean heat
    filename_ocheat = joinpath(@__DIR__, "data", "model_data", "sneasy_oceanheat_$(scenario)_$(Int(years[1]))_$(Int(years[end]))_$(Dates.format(now(),"dd-mm-yyyy")).csv")
    CSV.write(filename_ocheat, DataFrame([years,ocean_heat], ["Year","MAP Ocean Heat"]))

end

##==============================================================================
