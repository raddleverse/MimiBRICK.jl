using Mimi

#-------------------------------------------------------------------------------
# Create a function to run the DOECLIM+BRICK model over the historic period.
#-------------------------------------------------------------------------------

"""
    construct_run_doeclimbrick(calibration_start_year::Int, calibration_end_year::Int)

Create a function to run the DOECLIM+BRICK model over the historic period.
"""
function construct_run_doeclimbrick(calibration_start_year::Int, calibration_end_year::Int)

    # Load an instance of DOECLIM+BRICK model.
    # WARNING: for general use, use `m = create_brick_doeclim!(...[arguments here]...)`  instead
    m = Mimi.build(create_brick_doeclim(rcp_scenario="RCP85", start_year=calibration_start_year, end_year=calibration_end_year))

    # Get indices needed to normalize temperature anomalies relative to 1861-1880 mean (DOECLIM+BRICK starts in 1850 by default).
    temperature_norm_indices = findall((in)(1861:1880), 1850:calibration_end_year)

    # Get indices needed to normalize all sea level rise sources.
    sealevel_norm_indices_1961_1990 = findall((in)(1961:1990), 1850:calibration_end_year)
    sealevel_norm_indices_1992_2001 = findall((in)(1992:2001), 1850:calibration_end_year)

    # Given user settings, create a function to run DOECLIM+BRICK and return model output used for calibration.
    function run_doeclimbrick!(
        param::Array{Float64,1},
        modeled_temperature::Vector{Float64},
        modeled_ocean_heat::Vector{Float64},
        modeled_glaciers::Vector{Float64},
        modeled_greenland::Vector{Float64},
        modeled_antarctic::Vector{Float64},
        modeled_thermal_expansion::Vector{Float64},
        modeled_gmsl::Vector{Float64})

        # Assign names to uncertain model and initial condition parameters for convenience.
        # Note: This assumes "param" is the full vector of uncertain parameters with the same ordering as in "create_log_posterior_doeclim_brick.jl".
        temperature_0            = param[13]
        ocean_heat_0             = param[14]
        thermal_s???               = param[15]
        greenland_v???             = param[16]
        glaciers_v???              = param[17]
        glaciers_s???              = param[18]
        antarctic_s???             = param[19]
        heat_diffusivity         = param[20]
        rf_scale_aerosol         = param[21]
        ECS                      = param[22]
        thermal_??                = param[23]
        greenland_a              = param[24]
        greenland_b              = param[25]
        greenland_??              = param[26]
        greenland_??              = param[27]
        glaciers_?????              = param[28]
        glaciers_n               = param[29]
        anto_??                   = param[30]
        anto_??                   = param[31]
        antarctic_??              = param[32]
        antarctic_??              = param[33]
        antarctic_??              = param[34]
        antarctic_??              = param[35]
        antarctic_precip???        = param[36]
        antarctic_??              = param[37]
        antarctic_flow???          = param[38]
        antarctic_runoff_height??? = param[39]
        antarctic_c              = param[40]
        antarctic_bedheight???     = param[41]
        antarctic_slope          = param[42]
        antarctic_??              = param[43]
        antarctic_temp_threshold = param[44]

        #----------------------------------------------------------
        # Set DOECLIM+BRICK to use sampled parameter values.
        #----------------------------------------------------------

        # ---- Diffusion Ocean Energy balance CLIMate model (DOECLIM) ---- #
        update_param!(m, :doeclim, :t2co,  ECS)
        update_param!(m, :doeclim, :kappa, heat_diffusivity)

        # ---- Total Radiative Forcing ---- #
        update_param!(m, :radiativeforcing, :alpha, rf_scale_aerosol)

        # ----- Antarctic Ocean ----- #
        update_param!(m, :antarctic_ocean, :anto_??, anto_??)
        update_param!(m, :antarctic_ocean, :anto_??, anto_??)

        # ----- Antarctic Ice Sheet ----- #
        update_param!(m, :antarctic_icesheet, :ais_sea_level???,             antarctic_s???)
        update_param!(m, :antarctic_icesheet, :ais_bedheight???,             antarctic_bedheight???)
        update_param!(m, :antarctic_icesheet, :ais_slope,                  antarctic_slope)
        update_param!(m, :antarctic_icesheet, :ais_??,                      antarctic_??)
        update_param!(m, :antarctic_icesheet, :ais_runoffline_snowheight???, antarctic_runoff_height???)
        update_param!(m, :antarctic_icesheet, :ais_c,                      antarctic_c)
        update_param!(m, :antarctic_icesheet, :ais_precipitation???,         antarctic_precip???)
        update_param!(m, :antarctic_icesheet, :ais_??,                      antarctic_??)
        update_param!(m, :antarctic_icesheet, :ais_??,                      antarctic_??)
        update_param!(m, :antarctic_icesheet, :ais_iceflow???,               antarctic_flow???)
        update_param!(m, :antarctic_icesheet, :ais_??,                      antarctic_??)
        update_param!(m, :antarctic_icesheet, :ais_??,                      antarctic_??)
        update_param!(m, :antarctic_icesheet, :temperature_threshold,      antarctic_temp_threshold)
        update_param!(m, :antarctic_icesheet, :??,                          antarctic_??)

        # ----- Glaciers & Small Ice Caps ----- #
        update_param!(m, :glaciers_small_icecaps, :gsic_?????, glaciers_?????)
        update_param!(m, :glaciers_small_icecaps, :gsic_v???, glaciers_v???)
        update_param!(m, :glaciers_small_icecaps, :gsic_s???, glaciers_s???)
        update_param!(m, :glaciers_small_icecaps, :gsic_n,  glaciers_n)

        # ----- Greenland Ice Sheet ----- #
        update_param!(m, :greenland_icesheet, :greenland_a,  greenland_a)
        update_param!(m, :greenland_icesheet, :greenland_b,  greenland_b)
        update_param!(m, :greenland_icesheet, :greenland_??,  greenland_??)
        update_param!(m, :greenland_icesheet, :greenland_??,  greenland_??)
        update_param!(m, :greenland_icesheet, :greenland_v???, greenland_v???)

        # ----- Thermal Expansion ----- #
        update_param!(m, :thermal_expansion, :te_??,  thermal_??)
        update_param!(m, :thermal_expansion, :te_s???, thermal_s???)

        # Run model.
        run(m)

        #----------------------------------------------------------
        # Calculate model output being compared to observations.
        #----------------------------------------------------------

        # Global surface temperature anomaly (normalized to 1861-1880 mean with initial condition offset).
        modeled_temperature[:] = m[:doeclim, :temp] .- mean(m[:doeclim, :temp][temperature_norm_indices]) .+ temperature_0

        # Ocean heat content (with initial condition offset).
        modeled_ocean_heat[:] = m[:doeclim, :heat_mixed] .+ m[:doeclim, :heat_interior] .+ ocean_heat_0

        # Glaciers and small ice caps (normalized relative to 1961-1990 mean).
        modeled_glaciers[:] = m[:glaciers_small_icecaps, :gsic_sea_level] .- mean(m[:glaciers_small_icecaps, :gsic_sea_level][sealevel_norm_indices_1961_1990])

        # Greenland ice sheet (normalized relative to 1992-2001 ten year period to work with pooled data that includes IMBIE observations).
        modeled_greenland[:] = m[:greenland_icesheet, :greenland_sea_level] .- mean(m[:greenland_icesheet, :greenland_sea_level][sealevel_norm_indices_1961_1990])

        # Antarctic ice sheet (normalized relative to 1992-2001 ten year period to work with IMBIE data).
        modeled_antarctic[:] = m[:antarctic_icesheet, :ais_sea_level] .- mean(m[:antarctic_icesheet, :ais_sea_level][sealevel_norm_indices_1992_2001])

        # Sea level contribution from thermal expansion (calibrating to observed trends, so do not need to normalize).
        modeled_thermal_expansion[:] = m[:thermal_expansion, :te_sea_level]

		# Global mean sea level rise (normalize realtive to 1961-1990 mean).
		modeled_gmsl[:] = m[:global_sea_level, :sea_level_rise] .- mean(m[:global_sea_level, :sea_level_rise][sealevel_norm_indices_1961_1990])

        # Return results.
        return
    end

    # Return run model function.
    return run_doeclimbrick!
end

##------------------------------------------------------------------------------
## End
##------------------------------------------------------------------------------
