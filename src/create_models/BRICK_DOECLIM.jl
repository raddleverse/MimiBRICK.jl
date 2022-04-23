using CSVFiles
using DataFrames
using Distributions
using Mimi
using MimiSNEASY

# -------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------
# Function to create 'Building blocks for Relevant Ice and Climate Knowledge' (BRICK) model.
# -------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------

function create_brick_doeclim(;rcp_scenario::String = "RCP85", start_year::Int=1850, end_year::Int=2020)

    #-----------------------#
    # ----- Load Data ----- #
    #-----------------------#

    # Set model years.
	model_years = collect(start_year:end_year)

    # Find indices for BRICK start and end years relative to RCP time range of 1765-2500.
    index_start, index_end = findall((in)([start_year, end_year]), (1765:2500))

    # Load and clean up RCP radiative forcing data for DOEclim (options include "RCP26", "RCP45", "RCP60", and "RCP85").
    rcp_forcing_data = DataFrame(load(joinpath(@__DIR__, "..", "..", "data", "model_data", rcp_scenario*"_midyear_radforcings.csv"), skiplines_begin=58))[index_start:index_end, :]

    # Isolate radiative forcing from CO₂, aerosols, and all other sources.
    forcing_CO₂           = rcp_forcing_data.CO2_RF
    forcing_aerosols      = rcp_forcing_data.TOTAER_DIR_RF .+ rcp_forcing_data.CLOUD_TOT_RF
    forcing_other_sources = rcp_forcing_data.TOTAL_INCLVOLCANIC_RF .- forcing_CO₂ .- forcing_aerosols

    #-------------------------#
    # ----- Build BRICK ----- #
    #-------------------------#

    # Create a Mimi-model instance.
    brick_doeclim = Model()

    # Set model time horizon, defaults to 1850-2300.
    set_dimension!(brick_doeclim, :time, start_year:end_year)

    # Add in DOEclim and BRICK components.
    add_comp!(brick_doeclim, MimiSNEASY.radiativeforcing)
    add_comp!(brick_doeclim, MimiSNEASY.doeclim)
    add_comp!(brick_doeclim, antarctic_ocean)
    add_comp!(brick_doeclim, antarctic_icesheet)
    add_comp!(brick_doeclim, glaciers_small_icecaps)
    add_comp!(brick_doeclim, greenland_icesheet)
    add_comp!(brick_doeclim, thermal_expansion)
    add_comp!(brick_doeclim, landwater_storage)
    add_comp!(brick_doeclim, global_sea_level)

    #-------------------------------------#
    # ----- Assign Model Parameters ----- #
    #-------------------------------------#

    # ----- Shared Parameters ----- #

    add_shared_param!(brick_doeclim, :model_deltat, 1.0)
    connect_param!(brick_doeclim, :radiativeforcing, :deltat, :model_deltat)
    connect_param!(brick_doeclim, :doeclim, :deltat, :model_deltat)

    add_shared_param!(brick_doeclim, :model_seawater_freeze, -1.8)
    connect_param!(brick_doeclim, :antarctic_ocean, :seawater_freeze, :model_seawater_freeze)
    connect_param!(brick_doeclim, :antarctic_icesheet, :seawater_freeze, :model_seawater_freeze)

    # ----- Radiative Forcing ----- #

    update_param!(brick_doeclim, :radiativeforcing, :rf_co2, forcing_CO₂)
    update_param!(brick_doeclim, :radiativeforcing, :rf_aerosol, forcing_aerosols)
    update_param!(brick_doeclim, :radiativeforcing, :rf_other, forcing_other_sources)
    update_param!(brick_doeclim, :radiativeforcing, :alpha, 1.0)

    # ----- DOEclim ----- #

    update_param!(brick_doeclim, :doeclim, :t2co, 3.0)
    update_param!(brick_doeclim, :doeclim, :kappa, 1.0)

    # ----- Antarctic Ocean ----- #

    update_param!(brick_doeclim, :antarctic_ocean, :anto_α, 0.28)
    update_param!(brick_doeclim, :antarctic_ocean, :anto_β, 0.95)

    # ----- Antarctic Ice Sheet ----- #

    update_param!(brick_doeclim, :antarctic_icesheet, :ais_ρ_ice, 917.0)
    update_param!(brick_doeclim, :antarctic_icesheet, :ais_ρ_seawater, 1030.0)
    update_param!(brick_doeclim, :antarctic_icesheet, :ais_ρ_rock, 4000.0)
    update_param!(brick_doeclim, :antarctic_icesheet, :ais_sea_level₀, 0.0)
    update_param!(brick_doeclim, :antarctic_icesheet, :ais_ocean_temperature₀, 0.72)
    update_param!(brick_doeclim, :antarctic_icesheet, :ais_radius₀, 1.864e6)
    update_param!(brick_doeclim, :antarctic_icesheet, :ais_bedheight₀, 781.0)
    update_param!(brick_doeclim, :antarctic_icesheet, :ais_slope, 0.0006)
    update_param!(brick_doeclim, :antarctic_icesheet, :ais_μ, 11.0)
    update_param!(brick_doeclim, :antarctic_icesheet, :ais_runoffline_snowheight₀, 1400.0)
    update_param!(brick_doeclim, :antarctic_icesheet, :ais_c, 100.0)
    update_param!(brick_doeclim, :antarctic_icesheet, :ais_precipitation₀, 0.37)
    update_param!(brick_doeclim, :antarctic_icesheet, :ais_κ, 0.062)
    update_param!(brick_doeclim, :antarctic_icesheet, :ais_ν, 0.0086)
    update_param!(brick_doeclim, :antarctic_icesheet, :ais_iceflow₀, 1.2)
    update_param!(brick_doeclim, :antarctic_icesheet, :ais_γ, 2.9)
    update_param!(brick_doeclim, :antarctic_icesheet, :ais_α, 0.23)
    update_param!(brick_doeclim, :antarctic_icesheet, :ais_temperature_coefficient, 0.8365)
    update_param!(brick_doeclim, :antarctic_icesheet, :ais_temperature_intercept, 15.42)
    update_param!(brick_doeclim, :antarctic_icesheet, :ais_local_fingerprint, -1.18)
    update_param!(brick_doeclim, :antarctic_icesheet, :ocean_surface_area, 3.619e14)
    update_param!(brick_doeclim, :antarctic_icesheet, :temperature_threshold, -15.0)
    update_param!(brick_doeclim, :antarctic_icesheet, :λ, 0.0093)
    update_param!(brick_doeclim, :antarctic_icesheet, :include_ais_DSL, true)

    # ----- Glaciers & Small Ice Caps ----- #

    update_param!(brick_doeclim, :glaciers_small_icecaps, :gsic_β₀, 0.0013)
    update_param!(brick_doeclim, :glaciers_small_icecaps, :gsic_v₀, 0.376)
    update_param!(brick_doeclim, :glaciers_small_icecaps, :gsic_s₀, -0.0138)
    update_param!(brick_doeclim, :glaciers_small_icecaps, :gsic_n, 0.847)
    update_param!(brick_doeclim, :glaciers_small_icecaps, :gsic_teq, -0.15)

    # ----- Greenland Ice Sheet ----- #

    update_param!(brick_doeclim, :greenland_icesheet, :greenland_a, -1.37)
    update_param!(brick_doeclim, :greenland_icesheet, :greenland_b, 8.06)
    update_param!(brick_doeclim, :greenland_icesheet, :greenland_α, 0.0008)
    update_param!(brick_doeclim, :greenland_icesheet, :greenland_β, 0.00009)
    update_param!(brick_doeclim, :greenland_icesheet, :greenland_v₀, 7.52)

    # ----- Thermal Expansion ----- #

    update_param!(brick_doeclim, :thermal_expansion, :te_A, 3.619e14)
    update_param!(brick_doeclim, :thermal_expansion, :te_C, 3991.86795711963)
    update_param!(brick_doeclim, :thermal_expansion, :te_ρ, 1027.0)
    update_param!(brick_doeclim, :thermal_expansion, :te_α, 0.16)
    update_param!(brick_doeclim, :thermal_expansion, :te_s₀, 0.0)

    # ----- Landwater Storage ----- #

    update_param!(brick_doeclim, :landwater_storage, :lws₀, 0.0)
    update_param!(brick_doeclim, :landwater_storage, :first_projection_year, 2018)
    update_param!(brick_doeclim, :landwater_storage, :lws_random_sample, rand(Normal(0.0003, 0.00018), length(model_years)))


    #-----------------------------------------#
    #----- Create Component Connections ----- #
    #-----------------------------------------#

    # Set parameter connections (:component needing a value => :value name, :component producing the value => value name).

    # in default BRICK we use a model parameter :model_global_surface_temperature
    # here and connect all components to that, but here we will individually
    # connect the components to :doeclim since they are now pulling from another
    # component's variable and not from a model shared parameter
    connect_param!(brick_doeclim, :antarctic_icesheet => :global_surface_temperature,       :doeclim => :temp)
    connect_param!(brick_doeclim, :antarctic_ocean => :global_surface_temperature,          :doeclim => :temp)
    connect_param!(brick_doeclim, :glaciers_small_icecaps => :global_surface_temperature,   :doeclim => :temp)
    connect_param!(brick_doeclim, :greenland_icesheet => :global_surface_temperature,       :doeclim => :temp)

    connect_param!(brick_doeclim, :doeclim => :forcing, :radiativeforcing => :rf)

    connect_param!(brick_doeclim, :thermal_expansion => :ocean_heat_mixed,    :doeclim => :heat_mixed)
    connect_param!(brick_doeclim, :thermal_expansion => :ocean_heat_interior, :doeclim => :heat_mixed)

    connect_param!(brick_doeclim, :global_sea_level => :slr_glaciers_small_ice_caps, :glaciers_small_icecaps => :gsic_sea_level)
    connect_param!(brick_doeclim, :global_sea_level => :slr_greeland_icesheet,       :greenland_icesheet     => :greenland_sea_level)
    connect_param!(brick_doeclim, :global_sea_level => :slr_antartic_icesheet,       :antarctic_icesheet     => :ais_sea_level)
    connect_param!(brick_doeclim, :global_sea_level => :slr_thermal_expansion,       :thermal_expansion      => :te_sea_level)
    connect_param!(brick_doeclim, :global_sea_level => :slr_landwater_storage,       :landwater_storage      => :lws_sea_level)

    connect_param!(brick_doeclim, :antarctic_icesheet => :antarctic_ocean_temperature, :antarctic_ocean  => :anto_temperature)
    connect_param!(brick_doeclim, :antarctic_icesheet => :global_sea_level,            :global_sea_level => :sea_level_rise)

    # Return BRICK + DOEclim model.
    return brick_doeclim

end
