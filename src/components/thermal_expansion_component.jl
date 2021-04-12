@defcomp thermal_expansion begin

    # --------------------
    # Model Parameters
    # --------------------

    te_A                = Parameter()             # Global ocean surface area (m²).
    te_C                = Parameter()             # Heat capacity of conservative temperature (J kg⁻¹ °C⁻¹).
    te_ρ                = Parameter()             # Approximate density of global ocean (kg m⁻³).
    te_α                = Parameter()             # Global ocean-averaged thermal expansion coefficient (kg m⁻³ °C⁻¹).
    te_s₀               = Parameter()             # Initial sea level rise due to thermal expansion designated in 1850 (m SLE).
    ocean_heat_mixed    = Parameter(index=[time]) # Ocean heat content anomaly of mixed layer (10²² J).
    ocean_heat_interior = Parameter(index=[time]) # Ocean heat content anomaly of interior ocean (10²² J).

    # --------------------
    # Model Variables
    # --------------------

    ocean_heat_content = Variable(index=[time])   # Sum of ocean heat content in mixed layer and interior ocean.
    Δ_oceanheat        = Variable(index=[time])   # Change in total ocean heat content (J).
    te_sea_level       = Variable(index=[time])   # Cumulative sea level rise due to thermal expansion (m).

    # --------------------
    # Model Equations
    # --------------------

    function run_timestep(p, v, d, t)

        # Set initial conditions.
        if is_first(t)
            v.te_sea_level[t] = p.te_s₀
            v.Δ_oceanheat[t] = 0.0
            v.ocean_heat_content[t] = p.ocean_heat_mixed[t] + p.ocean_heat_interior[t]
        else
            # Calculate total ocean heat content anomaly.
            v.ocean_heat_content[t] = p.ocean_heat_mixed[t] + p.ocean_heat_interior[t]

            # Calculate change in total ocean heat content and convert to Joules.
            v.Δ_oceanheat[t] = (v.ocean_heat_content[t] - v.ocean_heat_content[t-1]) * 1e22

            # Calculate thermal expansion contribution to sea level rise given the change in global ocean heat content.
            v.te_sea_level[t] = v.te_sea_level[t-1] + v.Δ_oceanheat[t] * p.te_α / (p.te_A * p.te_C * p.te_ρ^2)
        end
    end
end
