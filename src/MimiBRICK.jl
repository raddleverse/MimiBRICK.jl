#using Distributions

module MimiBRICK

# Load required packages
using CSVFiles, DataFrames, Distributions, Mimi

# Load BRICK Mimi components.
include(joinpath("components", "antarctic_icesheet_component.jl"))
include(joinpath("components", "antarctic_ocean_component.jl"))
include(joinpath("components", "glaciers_small_icecaps_component.jl"))
include(joinpath("components", "global_sea_level_component.jl"))
include(joinpath("components", "greenland_icesheet_component.jl"))
include(joinpath("components", "landwater_storage_component.jl"))
include(joinpath("components", "thermal_expansion_component.jl"))

# -------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------
# Function to create 'Building blocks for Relevant Ice and Climate Knowledge' (BRICK) model.
# -------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------

function get_model(;rcp_scenario::String="RCP85", start_year::Int=1850, end_year::Int=2300)

    #-----------------------#
    # ----- Load Data ----- #
    #-----------------------#

    # Load exogenous time-series for global surface temperature and ocean heat content (output from SNEASY under RCP8.5).
    # NOTE: for now, only `rcp_scenario = "RCP85"` is supported
    temperature_scenario = DataFrame(load(joinpath(@__DIR__, "..", "data", "model_data", "sneasy_temperature_"*rcp_scenario*"_1850_2300.csv"), skiplines_begin=5))
    oceanheat_scenario   = DataFrame(load(joinpath(@__DIR__, "..", "data", "model_data", "sneasy_oceanheat_"*rcp_scenario*"_1850_2300.csv"), skiplines_begin=5))

    #-------------------------#
    # ----- Build BRICK ----- #
    #-------------------------#

    # Create a Mimi-model instance.
    brick = Model()

    # Set model time horizon, defaults to 1850-2300.
    set_dimension!(brick, :time, start_year:end_year)

    # Add in BRICK components.
    add_comp!(brick, antarctic_ocean)
    add_comp!(brick, antarctic_icesheet)
    add_comp!(brick, glaciers_small_icecaps)
    add_comp!(brick, greenland_icesheet)
    add_comp!(brick, thermal_expansion)
    add_comp!(brick, landwater_storage)
    add_comp!(brick, global_sea_level)

    #-------------------------------------#
    # ----- Assign Model Parameters ----- #
    #-------------------------------------#

    # ----- Antarctic Ocean ----- #

    set_param!(brick, :antarctic_ocean, :anto_α, 0.28)
    set_param!(brick, :antarctic_ocean, :anto_β, 0.95)

    # ----- Antarctic Ice Sheet ----- #

    set_param!(brick, :antarctic_icesheet, :ais_ρ_ice, 917.0)
    set_param!(brick, :antarctic_icesheet, :ais_ρ_seawater, 1030.0)
    set_param!(brick, :antarctic_icesheet, :ais_ρ_rock, 4000.0)
    set_param!(brick, :antarctic_icesheet, :ais_sea_level₀, 0.0)
    set_param!(brick, :antarctic_icesheet, :ais_ocean_temperature₀, 0.72)
    set_param!(brick, :antarctic_icesheet, :ais_radius₀, 1.864e6)
    set_param!(brick, :antarctic_icesheet, :ais_bedheight₀, 781.0)
    set_param!(brick, :antarctic_icesheet, :ais_slope, 0.0006)
    set_param!(brick, :antarctic_icesheet, :ais_μ, 11.0)
    set_param!(brick, :antarctic_icesheet, :ais_runoffline_snowheight₀, 1400.0)
    set_param!(brick, :antarctic_icesheet, :ais_c, 100.0)
    set_param!(brick, :antarctic_icesheet, :ais_precipitation₀, 0.37)
    set_param!(brick, :antarctic_icesheet, :ais_κ, 0.062)
    set_param!(brick, :antarctic_icesheet, :ais_ν, 0.0086)
    set_param!(brick, :antarctic_icesheet, :ais_iceflow₀, 1.2)
    set_param!(brick, :antarctic_icesheet, :ais_γ, 2.9)
    set_param!(brick, :antarctic_icesheet, :ais_α, 0.23)
    set_param!(brick, :antarctic_icesheet, :ais_temperature_coefficient, 0.8365)
    set_param!(brick, :antarctic_icesheet, :ais_temperature_intercept, 15.42)
    set_param!(brick, :antarctic_icesheet, :ais_local_fingerprint, -1.18)
    set_param!(brick, :antarctic_icesheet, :ocean_surface_area, 3.619e14)
    set_param!(brick, :antarctic_icesheet, :temperature_threshold, -15.0)
    set_param!(brick, :antarctic_icesheet, :λ, 0.0093)
    set_param!(brick, :antarctic_icesheet, :include_ais_DSL, true)

    # ----- Glaciers & Small Ice Caps ----- #

    set_param!(brick, :glaciers_small_icecaps, :gsic_β₀, 0.0013)
    set_param!(brick, :glaciers_small_icecaps, :gsic_v₀, 0.376)
    set_param!(brick, :glaciers_small_icecaps, :gsic_s₀, -0.0138)
    set_param!(brick, :glaciers_small_icecaps, :gsic_n, 0.847)
    set_param!(brick, :glaciers_small_icecaps, :gsic_teq, -0.15)

    # ----- Greenland Ice Sheet ----- #

    set_param!(brick, :greenland_icesheet, :greenland_a, -1.37)
    set_param!(brick, :greenland_icesheet, :greenland_b, 8.06)
    set_param!(brick, :greenland_icesheet, :greenland_α, 0.0008)
    set_param!(brick, :greenland_icesheet, :greenland_β, 0.00009)
    set_param!(brick, :greenland_icesheet, :greenland_v₀, 7.52)

    # ----- Thermal Expansion ----- #

    set_param!(brick, :thermal_expansion, :te_A, 3.619e14)
    set_param!(brick, :thermal_expansion, :te_C, 3991.86795711963)
    set_param!(brick, :thermal_expansion, :te_ρ, 1027.0)
    set_param!(brick, :thermal_expansion, :te_α, 0.16)
    set_param!(brick, :thermal_expansion, :te_s₀, 0.0)
    set_param!(brick, :thermal_expansion, :ocean_heat_mixed, zeros(length(start_year:end_year)))
    oceanheat_idx = findall((in)(start_year:end_year), oceanheat_scenario[!,:Year])
    set_param!(brick, :thermal_expansion, :ocean_heat_interior, oceanheat_scenario[oceanheat_idx, :Mean_Ocean_Heat])

    # ----- Landwater Storage ----- #

    set_param!(brick, :landwater_storage, :lws₀, 0.0)
    set_param!(brick, :landwater_storage, :first_projection_year, 2018)
    set_param!(brick, :landwater_storage, :lws_random_sample, rand(Normal(0.0003, 0.00018), length(start_year:end_year)))

    # ----- Set Parameters With Common Values Across Components ----- #

    set_param!(brick, :seawater_freeze, -1.8)
    temperature_idx = findall((in)(start_year:end_year), temperature_scenario[!,:Year])
    set_param!(brick, :global_surface_temperature, temperature_scenario[temperature_idx,:Mean_Temperature])

    #-----------------------------------------#
    #----- Create Component Connections ----- #
    #-----------------------------------------#

    # Set parameter connections (:component => :parameter).
    connect_param!(brick, :global_sea_level => :slr_glaciers_small_ice_caps, :glaciers_small_icecaps => :gsic_sea_level)
    connect_param!(brick, :global_sea_level => :slr_greeland_icesheet,       :greenland_icesheet     => :greenland_sea_level)
    connect_param!(brick, :global_sea_level => :slr_antartic_icesheet,       :antarctic_icesheet     => :ais_sea_level)
    connect_param!(brick, :global_sea_level => :slr_thermal_expansion,       :thermal_expansion      => :te_sea_level)
    connect_param!(brick, :global_sea_level => :slr_landwater_storage,       :landwater_storage      => :lws_sea_level)

    connect_param!(brick, :antarctic_icesheet => :antarctic_ocean_temperature, :antarctic_ocean  => :anto_temperature)
    connect_param!(brick, :antarctic_icesheet => :global_sea_level,            :global_sea_level => :sea_level_rise)

    # Return BRICK model.
    return brick
end

end # Module
