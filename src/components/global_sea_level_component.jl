using Mimi

@defcomp global_sea_level begin

    # --------------------
    # Model Parameters
    # --------------------

    slr_glaciers_small_ice_caps = Parameter(index=[time]) # Sea level contribution from glaciers and small ice caps (m).
    slr_greeland_icesheet       = Parameter(index=[time]) # Sea level contribution from Greenland ice sheet (m).
    slr_antartic_icesheet       = Parameter(index=[time]) # Sea level contribution from Antarctic ice sheet (m).
    slr_thermal_expansion       = Parameter(index=[time]) # Sea level contribution from thermal expansion (m).
    slr_landwater_storage       = Parameter(index=[time]) # Sea level contribution from landwater storage (m).

    # --------------------
    # Model Variables
    # --------------------

    sea_level_rise = Variable(index=[time])  # total sea level rise from all components (includes landwater storage for projection periods).

    # --------------------
    # Model Equations
    # --------------------

    function run_timestep(p, v, d, t)

        # Calculate global mean sea level rise as sum of sea level contributions from individual BRICK components.
        v.sea_level_rise[t] = p.slr_glaciers_small_ice_caps[t] + p.slr_greeland_icesheet[t] + p.slr_antartic_icesheet[t] + p.slr_thermal_expansion[t] + p.slr_landwater_storage[t]
    end
end
