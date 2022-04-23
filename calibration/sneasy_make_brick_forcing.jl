using CSV
using DataFrames
using MimiSNEASY
using Dates

##==============================================================================
## Read output from SNEASY RCP projections

dir_sneasybrick = joinpath(@__DIR__, "..", "results", "my_sneasybrick_results_20M_19-02-2022")

for rcp_scenario = ["RCP26","RCP45","RCP60","RCP85"]

    filename_sneasy_map = joinpath(dir_sneasybrick,"projections_csv",rcp_scenario,"projections_MAP_sneasybrick.csv")
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
    filename_temp = joinpath(@__DIR__, "..", "data", "model_data", "sneasy_temperature_$(rcp_scenario)_$(Int(years[1]))_$(Int(years[end]))_$(Dates.format(now(),"dd-mm-yyyy")).csv")
    CSV.write(filename_temp, DataFrame([years,temperature_norm], ["Year","MAP Temperature"]))
    # ocean heat
    filename_ocheat = joinpath(@__DIR__, "..", "data", "model_data", "sneasy_oceanheat_$(rcp_scenario)_$(Int(years[1]))_$(Int(years[end]))_$(Dates.format(now(),"dd-mm-yyyy")).csv")
    CSV.write(filename_ocheat, DataFrame([years,ocean_heat], ["Year","MAP Ocean Heat"]))

end

##==============================================================================
