using Mimi

#--------------------------------------------------------------------------------------------------------------
# Create a function to run the BRICK model over the historic period.
#--------------------------------------------------------------------------------------------------------------

function construct_run_brick(calibration_start_year::Int, calibration_end_year::Int)

    # Load an instance of DOECLIM+BRICK model.
    # WARNING: for general use, use `m = get_model!(...[arguments here]...)`  instead
    m = Mimi.build(get_model(rcp_scenario="RCP85", start_year=calibration_start_year, end_year=calibration_end_year))

    # Get indices needed to normalize temperature anomalies relative to 1861-1880 mean (DOECLIM+BRICK starts in 1850 by default).
    temperature_norm_indices = findall((in)(1861:1880), 1850:calibration_end_year)

    # Get indices needed to normalize all sea level rise sources.
    sealevel_norm_indices_1961_1990 = findall((in)(1961:1990), 1850:calibration_end_year)
    sealevel_norm_indices_1992_2001 = findall((in)(1992:2001), 1850:calibration_end_year)

    # Given user settings, create a function to run BRICK and return model output used for calibration.
    function run_brick!(
        param::Array{Float64,1},
        modeled_glaciers::Vector{Float64},
        modeled_greenland::Vector{Float64},
        modeled_antarctic::Vector{Float64},
        modeled_thermal_expansion::Vector{Float64},
        modeled_gmsl::Vector{Float64})

        # Assign names to uncertain model and initial condition parameters for convenience.
        # Note: This assumes "param" is the full vector of uncertain parameters with the same ordering as in "create_log_posterior_brick.jl".
        thermal_s₀               = param[9]
        greenland_v₀             = param[10]
        glaciers_v₀              = param[11]
        glaciers_s₀              = param[12]
        antarctic_s₀             = param[13]
        thermal_α                = param[14]
        greenland_a              = param[15]
        greenland_b              = param[16]
        greenland_α              = param[17]
        greenland_β              = param[18]
        glaciers_β₀              = param[19]
        glaciers_n               = param[20]
        anto_α                   = param[21]
        anto_β                   = param[22]
        antarctic_γ              = param[23]
        antarctic_α              = param[24]
        antarctic_μ              = param[25]
        antarctic_ν              = param[26]
        antarctic_precip₀        = param[27]
        antarctic_κ              = param[28]
        antarctic_flow₀          = param[29]
        antarctic_runoff_height₀ = param[30]
        antarctic_c              = param[31]
        antarctic_bedheight₀     = param[32]
        antarctic_slope          = param[33]
        antarctic_λ              = param[34]
        antarctic_temp_threshold = param[35]

        #----------------------------------------------------------
        # Set BRICK to use sampled parameter values.
        #----------------------------------------------------------

        # ----- Antarctic Ocean ----- #
        update_param!(m,  :antarctic_ocean, :anto_α, anto_α)
        update_param!(m,  :antarctic_ocean, :anto_β, anto_β)

        # ----- Antarctic Ice Sheet ----- #
        update_param!(m, :antarctic_icesheet, :ais_sea_level₀,             antarctic_s₀)
        update_param!(m, :antarctic_icesheet, :ais_bedheight₀,             antarctic_bedheight₀)
        update_param!(m, :antarctic_icesheet, :ais_slope,                  antarctic_slope)
        update_param!(m, :antarctic_icesheet, :ais_μ,                      antarctic_μ)
        update_param!(m, :antarctic_icesheet, :ais_runoffline_snowheight₀, antarctic_runoff_height₀)
        update_param!(m, :antarctic_icesheet, :ais_c,                      antarctic_c)
        update_param!(m, :antarctic_icesheet, :ais_precipitation₀,         antarctic_precip₀)
        update_param!(m, :antarctic_icesheet, :ais_κ,                      antarctic_κ)
        update_param!(m, :antarctic_icesheet, :ais_ν,                      antarctic_ν)
        update_param!(m, :antarctic_icesheet, :ais_iceflow₀,               antarctic_flow₀)
        update_param!(m, :antarctic_icesheet, :ais_γ,                      antarctic_γ)
        update_param!(m, :antarctic_icesheet, :ais_α,                      antarctic_α)
        update_param!(m, :antarctic_icesheet, :temperature_threshold,      antarctic_temp_threshold)
        update_param!(m, :antarctic_icesheet, :λ,                          antarctic_λ)

        # ----- Glaciers & Small Ice Caps ----- #
        update_param!(m, :glaciers_small_icecaps, :gsic_β₀, glaciers_β₀)
        update_param!(m, :glaciers_small_icecaps, :gsic_v₀, glaciers_v₀)
        update_param!(m, :glaciers_small_icecaps, :gsic_s₀, glaciers_s₀)
        update_param!(m, :glaciers_small_icecaps, :gsic_n,  glaciers_n)

        # ----- Greenland Ice Sheet ----- #
        update_param!(m, :greenland_icesheet, :greenland_a,  greenland_a)
        update_param!(m, :greenland_icesheet, :greenland_b,  greenland_b)
        update_param!(m, :greenland_icesheet, :greenland_α,  greenland_α)
        update_param!(m, :greenland_icesheet, :greenland_β,  greenland_β)
        update_param!(m, :greenland_icesheet, :greenland_v₀, greenland_v₀)

        # ----- Thermal Expansion ----- #
        update_param!(m, :thermal_expansion, :te_α,  thermal_α)
        update_param!(m, :thermal_expansion, :te_s₀, thermal_s₀)

        # Run model.
        run(m)

        #----------------------------------------------------------
        # Calculate model output being compared to observations.
        #----------------------------------------------------------

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
    return run_brick!
end
