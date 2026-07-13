using Mimi

#-------------------------------------------------------------------------------
# Create a function to run the SNEASY+BRICK model over the historic period.
#-------------------------------------------------------------------------------

"""
    construct_run_sneasybrick(calibration_start_year::Int, calibration_end_year::Int, glacier_model::Symbol)

Create a function to run the SNEASY+BRICK model over the historic period.
"""
function construct_run_sneasybrick(calibration_start_year::Int, calibration_end_year::Int; glacier_model::Symbol=:mengel)

    # Load an instance of SNEASY+BRICK model.
    # WARNING: for general use, use `m = create_sneasy_brick!(...[arguments here]...)`  instead
    m = Mimi.build(create_sneasy_brick(; ssprcp_scenario = "ssp245", start_year=calibration_start_year, end_year=calibration_end_year, glacier_model=glacier_model))

    # Get indices needed to normalize temperature anomalies relative to 1861-1880 mean (SNEASY+BRICK starts in 1850 by default).
    temperature_norm_indices = findall((in)(1861:1880), 1850:calibration_end_year)

    # Get indices needed to normalize all sea level rise sources.
    sealevel_norm_indices_1961_1990 = findall((in)(1961:1990), 1850:calibration_end_year)
    sealevel_norm_indices_1992_2001 = findall((in)(1992:2001), 1850:calibration_end_year)

    # Given user settings, create a function to run SNEASY+BRICK and return model output used for calibration.
    function run_sneasybrick!(
        param::Array{Float64,1},
        modeled_CO₂::Vector{Float64},
        modeled_oceanCO₂_flux::Vector{Float64},
        modeled_temperature::Vector{Float64},
        modeled_ocean_heat::Vector{Float64},
        modeled_glaciers::Vector{Float64},
        modeled_greenland::Vector{Float64},
        modeled_antarctic::Vector{Float64},
        modeled_thermal_expansion::Vector{Float64},
        modeled_gmsl::Vector{Float64})

        # Assign names to uncertain model and initial condition parameters for convenience.
        # Note: This assumes "param" is the full vector of uncertain parameters with the same ordering as in "create_log_posterior_sneasy_brick.jl".
        if glacier_model==:gsic
            CO₂_0                    = param[15]
            N₂O_0                    = param[16]
            temperature_0            = param[17]
            ocean_heat_0             = param[18]
            thermal_s₀               = param[19]
            greenland_v₀             = param[20]
            glaciers_v₀              = param[21]
            glaciers_s₀              = param[22]
            antarctic_s₀             = param[23]
            Q10                      = param[24]
            CO₂_fertilization        = param[25]
            CO₂_diffusivity          = param[26]
            heat_diffusivity         = param[27]
            rf_scale_aerosol         = param[28]
            ECS                      = param[29]
            thermal_α                = param[30]
            greenland_a              = param[31]
            greenland_b              = param[32]
            greenland_α              = param[33]
            greenland_β              = param[34]
            glaciers_β₀              = param[35]
            glaciers_n               = param[36]
            anto_α                   = param[37]
            anto_β                   = param[38]
            antarctic_γ              = param[39]
            antarctic_α              = param[40]
            antarctic_μ              = param[41]
            antarctic_ν              = param[42]
            antarctic_precip₀        = param[43]
            antarctic_κ              = param[44]
            antarctic_flow₀          = param[45]
            antarctic_runoff_height₀ = param[46]
            antarctic_c              = param[47]
            antarctic_bedheight₀     = param[48]
            antarctic_slope          = param[49]
            antarctic_λ              = param[50]
            antarctic_temp_threshold = param[51]
        elseif glacier_model==:mengel
            CO₂_0                    = param[15]
            N₂O_0                    = param[16]
            temperature_0            = param[17]
            ocean_heat_0             = param[18]
            thermal_s₀               = param[19]
            greenland_v₀             = param[20]
            glaciers_sl0             = param[21]
            antarctic_s₀             = param[22]
            Q10                      = param[23]
            CO₂_fertilization        = param[24]
            CO₂_diffusivity          = param[25]
            heat_diffusivity         = param[26]
            rf_scale_aerosol         = param[27]
            ECS                      = param[28]
            thermal_α                = param[29]
            greenland_a              = param[30]
            greenland_b              = param[31]
            greenland_α              = param[32]
            greenland_β              = param[33]
            glaciers_a               = param[34] # update
            glaciers_b               = param[35]
            glaciers_T_lia           = param[36]
            glaciers_f               = param[37]
            glaciers_tau_fast        = param[38]
            glaciers_tau_slow        = param[39]
            anto_α                   = param[40]
            anto_β                   = param[41]
            antarctic_γ              = param[42]
            antarctic_α              = param[43]
            antarctic_μ              = param[44]
            antarctic_ν              = param[45]
            antarctic_precip₀        = param[46]
            antarctic_κ              = param[47]
            antarctic_flow₀          = param[48]
            antarctic_runoff_height₀ = param[49]
            antarctic_c              = param[50]
            antarctic_bedheight₀     = param[51]
            antarctic_slope          = param[52]
            antarctic_λ              = param[53]
            antarctic_temp_threshold = param[54]
        end

        #----------------------------------------------------------
        # Set SNEASY+BRICK to use sampled parameter values.
        #----------------------------------------------------------

        # ---- Shared Parameters ---- #
        update_param!(m, :model_CO₂_0, CO₂_0) # connects to :ccm and :rfco2

        # ---- Diffusion Ocean Energy balance CLIMate model (DOECLIM) ---- #
        update_param!(m, :doeclim, :t2co,  ECS)
        update_param!(m, :doeclim, :kappa, heat_diffusivity)

        # ---- Carbon Cycle ---- #
        update_param!(m, :ccm, :Q10,     Q10)
        update_param!(m, :ccm, :Beta,    CO₂_fertilization)
        update_param!(m, :ccm, :Eta,     CO₂_diffusivity)

        # ---- Carbon Dioxide Radiative Forcing ---- #
        update_param!(m, :rfco2, :N₂O_0, N₂O_0)

        # ---- Total Radiative Forcing ---- #
        update_param!(m, :radiativeforcing, :alpha, rf_scale_aerosol)

        # ----- Antarctic Ocean ----- #
        update_param!(m, :antarctic_ocean, :anto_α, anto_α)
        update_param!(m, :antarctic_ocean, :anto_β, anto_β)

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
        if glacier_model==:gsic
            update_param!(m, :glaciers_small_icecaps, :gsic_β₀, glaciers_β₀)
            update_param!(m, :glaciers_small_icecaps, :gsic_v₀, glaciers_v₀)
            update_param!(m, :glaciers_small_icecaps, :gsic_s₀, glaciers_s₀)
            update_param!(m, :glaciers_small_icecaps, :gsic_n,  glaciers_n)
        elseif glacier_model==:mengel
            update_param!(m, :glaciers_small_icecaps, :gic_a, glaciers_a)
            update_param!(m, :glaciers_small_icecaps, :gic_b, glaciers_b)
            update_param!(m, :glaciers_small_icecaps, :gic_T_lia, glaciers_T_lia)
            update_param!(m, :glaciers_small_icecaps, :gic_f, glaciers_f)
            update_param!(m, :glaciers_small_icecaps, :gic_tau_fast, glaciers_tau_fast)
            update_param!(m, :glaciers_small_icecaps, :gic_tau_slow, glaciers_tau_slow)
            update_param!(m, :glaciers_small_icecaps, :gic_sl0, glaciers_sl0)
        end

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

        # Atmospheric concentration of CO₂.
        #modeled_CO₂[:] = m[:ccm, :CO₂_0]
        modeled_CO₂[:] = m[:ccm, :atmco2]

        # Global surface temperature anomaly (normalized to 1861-1880 mean with initial condition offset).
        modeled_temperature[:] = m[:doeclim, :temp] .- mean(m[:doeclim, :temp][temperature_norm_indices]) .+ temperature_0

        # Ocean carbon flux (Note: timesteps cause last `atm_oc_flux` value to equal `missing`, so exclude it here).
        modeled_oceanCO₂_flux[1:end-1] = m[:ccm, :atm_oc_flux][1:end-1]

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
    return run_sneasybrick!
end

##------------------------------------------------------------------------------
## End
##------------------------------------------------------------------------------
