#--------------------------------------------------------------------------------------------------------------
# Create a function to run the BRICK model over the historic period.
#--------------------------------------------------------------------------------------------------------------

# Load model file.
using MimiBRICK

function construct_run_brick(calibration_start_year::Int, calibration_end_year::Int)

    # Load an instance of DOECLIM+BRICK model.
    #m = MimiBRICK_DOECLIM.create_brick_doeclim(rcp_scenario="RCP85", start_year=calibration_start_year, end_year=calibration_end_year)
    m = MimiBRICK.get_model() # TODO - add optional arguments for RCP scenario, beg/end years?

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
        #temperature_0            = param[13]
        #ocean_heat_0             = param[14]
        thermal_s₀               = param[15]
        greenland_v₀             = param[16]
        glaciers_v₀              = param[17]
        glaciers_s₀              = param[18]
        antarctic_s₀             = param[19]
        #heat_diffusivity         = param[20]
        #rf_scale_aerosol         = param[21]
        #ECS                      = param[22]
        thermal_α                = param[23]
        greenland_a              = param[24]
        greenland_b              = param[25]
        greenland_α              = param[26]
        greenland_β              = param[27]
        glaciers_β₀              = param[28]
        glaciers_n               = param[29]
        anto_α                   = param[30]
        anto_β                   = param[31]
        antarctic_γ              = param[32]
        antarctic_α              = param[33]
        antarctic_μ              = param[34]
        antarctic_ν              = param[35]
        antarctic_precip₀        = param[36]
        antarctic_κ              = param[37]
        antarctic_flow₀          = param[38]
        antarctic_runoff_height₀ = param[39]
        antarctic_c              = param[40]
        antarctic_bedheight₀     = param[41]
        antarctic_slope          = param[42]
        antarctic_λ              = param[43]
        antarctic_temp_threshold = param[44]

        #----------------------------------------------------------
        # Set BRICK to use sampled parameter values.
        #----------------------------------------------------------

        # ----- Antarctic Ocean ----- #
        update_param!(m, :anto_α, anto_α)
        update_param!(m, :anto_β, anto_β)

        # ----- Antarctic Ice Sheet ----- #
        update_param!(m, :ais_sea_level₀,             antarctic_s₀)
        update_param!(m, :ais_bedheight₀,             antarctic_bedheight₀)
        update_param!(m, :ais_slope,                  antarctic_slope)
        update_param!(m, :ais_μ,                      antarctic_μ)
        update_param!(m, :ais_runoffline_snowheight₀, antarctic_runoff_height₀)
        update_param!(m, :ais_c,                      antarctic_c)
        update_param!(m, :ais_precipitation₀,         antarctic_precip₀)
        update_param!(m, :ais_κ,                      antarctic_κ)
        update_param!(m, :ais_ν,                      antarctic_ν)
        update_param!(m, :ais_iceflow₀,               antarctic_flow₀)
        update_param!(m, :ais_γ,                      antarctic_γ)
        update_param!(m, :ais_α,                      antarctic_α)
        update_param!(m, :temperature_threshold,      antarctic_temp_threshold)
        update_param!(m, :λ,                          antarctic_λ)

        # ----- Glaciers & Small Ice Caps ----- #
        update_param!(m, :gsic_β₀, glaciers_β₀)
        update_param!(m, :gsic_v₀, glaciers_v₀)
        update_param!(m, :gsic_s₀, glaciers_s₀)
        update_param!(m, :gsic_n,  glaciers_n)

        # ----- Greenland Ice Sheet ----- #
        update_param!(m, :greenland_a,  greenland_a)
        update_param!(m, :greenland_b,  greenland_b)
        update_param!(m, :greenland_α,  greenland_α)
        update_param!(m, :greenland_β,  greenland_β)
        update_param!(m, :greenland_v₀, greenland_v₀)

        # ----- Thermal Expansion ----- #
        update_param!(m, :te_α,  thermal_α)
        update_param!(m, :te_s₀, thermal_s₀)

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
