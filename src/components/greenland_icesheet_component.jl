@defcomp greenland_icesheet begin

    # --------------------
    # Model Parameters
    # --------------------

    greenland_a                 = Parameter()             # Sensitivity of the equilibrium volume to changes in temperature (mSLE °C⁻¹).
    greenland_b                 = Parameter()             # Equilibrium volume for 0 °C global surface temperature anomaly (mSLE).
    greenland_α                 = Parameter()             # Temperature sensitivity of Greenland ice sheet exponential decay rate (yr⁻¹ °C⁻¹).
    greenland_β                 = Parameter()             # Exponential decay rate for 0 °C global surface temperature anomaly (yr⁻¹).
    greenland_v₀                = Parameter()             # Initial volume of Greenland ice sheet (m SLE).
    global_surface_temperature  = Parameter(index=[time]) # Global mean surface temperature relative to pre-industrial (°C).

    # --------------------
    # Model Variables
    # --------------------
    τ_inv               = Variable(index=[time]) # E-folding timescale of Greenland ice sheet volume changes due to changes in global temperature (yr⁻¹).
    eq_volume           = Variable(index=[time]) # Equilibrium ice sheet volume where sea level contribution is 0 (mSLE).
    greenland_volume    = Variable(index=[time]) # Volume of Greenland ice sheet (mSLE)
    greenland_sea_level = Variable(index=[time]) # Cumulative sea level rise contribution from Greenland ice sheet (m).

    # --------------------
    # Model Equations
    # --------------------

    function run_timestep(p, v, d, t)

        # Calculate equilibrium ice sheet volume for current temperature anomaly.
        v.eq_volume[t] = p.greenland_a * p.global_surface_temperature[t] + p.greenland_b

        if is_first(t)

            # Set initial conditions for volume.
            v.greenland_volume[t] = p.greenland_v₀

            # Calculate e-folding time scale based on initial conditions.
            v.τ_inv[t] = p.greenland_α * p.global_surface_temperature[t] + p.greenland_β

        else

            # Calculate current Greenland ice sheet volume.
            v.greenland_volume[t] = v.greenland_volume[t-1] + v.τ_inv[t-1] * (v.eq_volume[t-1] - v.greenland_volume[t-1])

            # Calculate new e-folding time scale dependent on temperature and remaining ice sheet volume.
            v.τ_inv[t] = (p.greenland_α * p.global_surface_temperature[t] + p.greenland_β) * (v.greenland_volume[t] / p.greenland_v₀)

        end

        # Calculate total sea level contribution from Greenland Ice Sheet (with a check in case the ice sheet fully melts).
        if v.greenland_volume[t] > 0.0
            v.greenland_sea_level[t] = p.greenland_v₀- v.greenland_volume[t]
        else
            v.greenland_sea_level[t] = p.greenland_v₀
        end
    end
end
