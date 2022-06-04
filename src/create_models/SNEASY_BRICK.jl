using CSVFiles
using DataFrames
using Distributions
using Mimi
using MimiSNEASY

#-------------------------------------------------------------------------------
# Function to run SNEASY-BRICK climate model over historic period.
# ------------------------------------------------------------------------------
"""
    create_sneasy_brick(;rcp_scenario::String = "RCP85", start_year::Int=1850, end_year::Int=2020)

RETURN A MIMI MODEL INSTANCE WITH MIMIBRICK AND MIMISNEASY COUPLED TOGETHER

Description: This function loads forcing data, sets up model parameters, and
makes the model component variable connections.

Function Arguments:

        rcp_scenario = RCP scenario for exogenous forcing
        start_year   = initial year of the simulation period
        end_year     = ending year of the simulation period
"""
function create_sneasy_brick(; rcp_scenario::String="RCP85", start_year::Int=1850, end_year::Int=2020)

 	# ---------------------------------------------
    # Load and clean up necessary data.
    # ---------------------------------------------

	# Set model years.
	model_years = collect(start_year:end_year)

    # Find indices for RCP data (1765-2500) corresponding to DICE years.
    rcp_indices = findall((in)(model_years), 1765:2500)

	# Load emissions and forcing data (index into appropriate years).
  	rcp_emissions      = DataFrame(load(joinpath(@__DIR__, "..", "..", "data", "model_data", rcp_scenario*"_emissions.csv"), skiplines_begin=36))
    rcp_concentrations = DataFrame(load(joinpath(@__DIR__, "..", "..", "data", "model_data", rcp_scenario*"_concentrations.csv"), skiplines_begin=37))
    rcp_forcing        = DataFrame(load(joinpath(@__DIR__, "..", "..", "data", "model_data", rcp_scenario*"_midyear_radforcings.csv"), skiplines_begin=58))

    # Calculate CO₂ emissions.
    rcp_co2_emissions = (rcp_emissions.FossilCO2 .+ rcp_emissions.OtherCO2)[rcp_indices]

    # Get RCP N₂O concentrations (used for CO₂ radiative forcing calculations).
    rcp_n2o_concentration = rcp_concentrations.N2O[rcp_indices]

    # Calculate exogenous radiative forcings.
    rcp_aerosol_forcing    = (rcp_forcing.TOTAER_DIR_RF .+ rcp_forcing.CLOUD_TOT_RF)[rcp_indices]
    rcp_other_forcing      = (rcp_forcing.TOTAL_INCLVOLCANIC_RF .- rcp_forcing.CO2_RF .- rcp_forcing.TOTAER_DIR_RF .- rcp_forcing.CLOUD_TOT_RF)[rcp_indices]

  	# Get an instance of Mimi-BRICK sea level rise model and set time dimension.
	m = MimiSNEASY.get_model(start_year=start_year, end_year=end_year)

	# Add in BRICK sea level rise components.
	add_comp!(m, antarctic_ocean,           after = :ccm)
	add_comp!(m, antarctic_icesheet,        after = :antarctic_ocean)
	add_comp!(m, glaciers_small_icecaps,    after = :antarctic_icesheet)
	add_comp!(m, greenland_icesheet,        after = :glaciers_small_icecaps)
	add_comp!(m, thermal_expansion,         after = :greenland_icesheet)
	add_comp!(m, landwater_storage,         after = :thermal_expansion)
	add_comp!(m, global_sea_level,          after = :landwater_storage)

    # Set new time dimension for all coupled model components.
    set_dimension!(m, :time, model_years)

	# Set all BRICK parameters and create model connections.

	# ----- Antarctic Ocean ----- #

    update_param!(m, :antarctic_ocean, :anto_α, 0.28)
    update_param!(m, :antarctic_ocean, :anto_β, 0.95)

    # ----- Antarctic Ice Sheet ----- #

    update_param!(m, :antarctic_icesheet, :ais_ρ_ice, 917.0)
    update_param!(m, :antarctic_icesheet, :ais_ρ_seawater, 1030.0)
    update_param!(m, :antarctic_icesheet, :ais_ρ_rock, 4000.0)
    update_param!(m, :antarctic_icesheet, :ais_sea_level₀, 0.0)
    update_param!(m, :antarctic_icesheet, :ais_ocean_temperature₀, 0.72)
    update_param!(m, :antarctic_icesheet, :ais_radius₀, 1.864e6)
    update_param!(m, :antarctic_icesheet, :ais_bedheight₀, 781.0)
    update_param!(m, :antarctic_icesheet, :ais_slope, 0.0006)
    update_param!(m, :antarctic_icesheet, :ais_μ, 11.0)
    update_param!(m, :antarctic_icesheet, :ais_runoffline_snowheight₀, 1400.0)
    update_param!(m, :antarctic_icesheet, :ais_c, 100.0)
    update_param!(m, :antarctic_icesheet, :ais_precipitation₀, 0.37)
    update_param!(m, :antarctic_icesheet, :ais_κ, 0.062)
    update_param!(m, :antarctic_icesheet, :ais_ν, 0.0086)
    update_param!(m, :antarctic_icesheet, :ais_iceflow₀, 1.2)
    update_param!(m, :antarctic_icesheet, :ais_γ, 2.9)
    update_param!(m, :antarctic_icesheet, :ais_α, 0.23)
    update_param!(m, :antarctic_icesheet, :ais_temperature_coefficient, 0.8365)
    update_param!(m, :antarctic_icesheet, :ais_temperature_intercept, 15.42)
    update_param!(m, :antarctic_icesheet, :ais_local_fingerprint, -1.18)
    update_param!(m, :antarctic_icesheet, :ocean_surface_area, 3.619e14)
    update_param!(m, :antarctic_icesheet, :temperature_threshold, -15.0)
    update_param!(m, :antarctic_icesheet, :λ, 0.0093)
    update_param!(m, :antarctic_icesheet, :include_ais_DSL, true)

    # ----- Glaciers & Small Ice Caps ----- #

    update_param!(m, :glaciers_small_icecaps, :gsic_β₀, 0.0013)
    update_param!(m, :glaciers_small_icecaps, :gsic_v₀, 0.376)
    update_param!(m, :glaciers_small_icecaps, :gsic_s₀, -0.0138)
    update_param!(m, :glaciers_small_icecaps, :gsic_n, 0.847)
    update_param!(m, :glaciers_small_icecaps, :gsic_teq, -0.15)

    # ----- Greenland Ice Sheet ----- #

    update_param!(m, :greenland_icesheet, :greenland_a, -1.37)
    update_param!(m, :greenland_icesheet, :greenland_b, 8.06)
    update_param!(m, :greenland_icesheet, :greenland_α, 0.0008)
    update_param!(m, :greenland_icesheet, :greenland_β, 0.00009)
    update_param!(m, :greenland_icesheet, :greenland_v₀, 7.52)

    # ----- Thermal Expansion ----- #

    update_param!(m, :thermal_expansion, :te_A, 3.619e14)
    update_param!(m, :thermal_expansion, :te_C, 3991.86795711963)
    update_param!(m, :thermal_expansion, :te_ρ, 1027.0)
    update_param!(m, :thermal_expansion, :te_α, 0.16)
    update_param!(m, :thermal_expansion, :te_s₀, 0.0)

    # ----- Landwater Storage ----- #

    update_param!(m, :landwater_storage, :lws₀, 0.0)
    update_param!(m, :landwater_storage, :first_projection_year, 2018)
    update_param!(m, :landwater_storage, :lws_random_sample, rand(Normal(0.0003, 0.00018), length(model_years)))

    # ----- Set Parameters With Common Values Across Components ----- #

    add_shared_param!(m, :model_seawater_freeze, -1.8)
    connect_param!(m, :antarctic_icesheet, :seawater_freeze, :model_seawater_freeze)
    connect_param!(m, :antarctic_ocean, :seawater_freeze, :model_seawater_freeze)

    # ----- SNEASY RCP Scenario Specific Parameters ----- #

	update_param!(m, :ccm, :CO2_emissions, rcp_co2_emissions)
	update_param!(m, :rfco2, :N₂O, rcp_n2o_concentration)
 	update_param!(m, :radiativeforcing, :rf_aerosol, rcp_aerosol_forcing)
 	update_param!(m, :radiativeforcing, :rf_other, rcp_other_forcing)

    #-----------------------------------------#
    #----- Create Component Connections ----- #
    #-----------------------------------------#

    # Set parameter connections (:component => :parameter).

    # in default BRICK we use a model parameter :model_global_surface_temperature
    # here and connect all components to that, but here we will individually
    # connect the components to :doeclim since they are now pulling from another
    # component's variable and not from a model shared parameter
    connect_param!(m, :antarctic_icesheet =>        :global_surface_temperature, :doeclim => :temp)
    connect_param!(m, :antarctic_ocean =>           :global_surface_temperature, :doeclim => :temp)
    connect_param!(m, :glaciers_small_icecaps =>    :global_surface_temperature, :doeclim => :temp)
    connect_param!(m, :greenland_icesheet =>        :global_surface_temperature, :doeclim => :temp)

    connect_param!(m, :thermal_expansion => :ocean_heat_mixed,    :doeclim => :heat_mixed)
    connect_param!(m, :thermal_expansion => :ocean_heat_interior, :doeclim => :heat_interior)

    connect_param!(m, :global_sea_level => :slr_glaciers_small_ice_caps, :glaciers_small_icecaps => :gsic_sea_level)
    connect_param!(m, :global_sea_level => :slr_greeland_icesheet,       :greenland_icesheet     => :greenland_sea_level)
    connect_param!(m, :global_sea_level => :slr_antartic_icesheet,       :antarctic_icesheet     => :ais_sea_level)
    connect_param!(m, :global_sea_level => :slr_thermal_expansion,       :thermal_expansion      => :te_sea_level)
    connect_param!(m, :global_sea_level => :slr_landwater_storage,       :landwater_storage      => :lws_sea_level)

    connect_param!(m, :antarctic_icesheet => :antarctic_ocean_temperature, :antarctic_ocean  => :anto_temperature)
    connect_param!(m, :antarctic_icesheet => :global_sea_level,            :global_sea_level => :sea_level_rise)

    # Return SNEASY-BRICK model.
    return m

end
