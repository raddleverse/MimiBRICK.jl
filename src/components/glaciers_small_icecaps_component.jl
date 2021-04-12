@defcomp glaciers_small_icecaps begin

    # --------------------
    # Model Parameters
    # --------------------

    gsic_β₀                    = Parameter()             # Initial mass-balance sensitivity to global surface temperature (m °C⁻¹ yr⁻¹).
    gsic_v₀                    = Parameter()             # Initial total volume of glaciers and small ice caps (m SLE).
    gsic_s₀                    = Parameter()             # Initial cumulative sea level rise from glaciers and small ice caps (m).
    gsic_n                     = Parameter()             # Exponential parameter for area to volume scaling.
    gsic_teq                   = Parameter()             # Theoretical equilibrium temperture at which glacier and small ice cap mass balance is in steady state (°C).
    global_surface_temperature = Parameter(index=[time]) # Global mean surface temperature anomaly relative to pre-industrial (°C).

    # --------------------
    # Model Variables
    # --------------------

    gsic_sea_level = Variable(index=[time])   # Cumulative sea level rise from glaciers and small ice caps (m).

    # --------------------
    # Model Equations
    # --------------------

    function run_timestep(p, v, d, t)

        if is_first(t)

            # Set initial sea level contribution from glaciers and small ice caps.
            v.gsic_sea_level[t] = p.gsic_s₀

        else

            # Calculate sea level contribution from glaciers and small ice caps (with a check in case they have fully melted).
            if v.gsic_sea_level[t-1] < p.gsic_v₀
                v.gsic_sea_level[t] = v.gsic_sea_level[t-1] + p.gsic_β₀ * (p.global_surface_temperature[t-1] - p.gsic_teq)*(1 - v.gsic_sea_level[t-1]/p.gsic_v₀)^p.gsic_n
            else
                v.gsic_sea_level[t] = p.gsic_v₀
            end
        end
    end
end
