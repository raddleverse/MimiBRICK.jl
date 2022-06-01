using Mimi

@defcomp antarctic_ocean begin

    # --------------------
    # Model Parameters
    # --------------------

    anto_α                     = Parameter()             # Sensitivity of the Antarctic Ocean temperature to global surface temperature (°C °C⁻¹).
    anto_β                     = Parameter()             # Antarctic ocean temperature for 0 °C global surface temperature anomaly (°C).
    seawater_freeze            = Parameter()             # Freezing temperature of seawater (°C).
    global_surface_temperature = Parameter(index=[time]) # Global mean surface temperature relative to pre-industrial (°C).

    # --------------------
    # Model Variables
    # --------------------

    anto_temperature            = Variable(index=[time]) # Ocean surface temperature at Antarctica (°C).

    # --------------------
    # Model Equations
    # --------------------

    function run_timestep(p, v, d, t)

        # Define temporary numerator and denominator terms (just for convenience).
        term_1 = p.anto_α * p.global_surface_temperature[t] + p.anto_β - p.seawater_freeze
        term_2 = 1 + exp( - (p.anto_α * p.global_surface_temperature[t] + p.anto_β - p.seawater_freeze) / p.anto_α)

        # Calculate antarctic ocean temperature (bounded below by seawater freezing point).
        v.anto_temperature[t] = p.seawater_freeze + (term_1 / term_2)
    end
end
