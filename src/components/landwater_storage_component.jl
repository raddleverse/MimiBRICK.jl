@defcomp landwater_storage begin

    # --------------------
    # Model Parameters
    # --------------------

    lws₀                  = Parameter()             # Initial landwater storage trend value applying to first projection year (m).
    first_projection_year = Parameter{Int}()        # First model projection year (necessary because calibrating to Church & White data that removed landwater storage trend).
    lws_distribution      = Parameter(index=[time]) # Samples from Normal distribution for annual land water storage trend.
    #lws_distribution     = Parameter{Distributions.Normal{Float64}}() # Normal distribution to sample for annual land water storage trend.

    # --------------------
    # Model Variables
    # --------------------

    lws_sea_level = Variable(index=[time])  # Cumulative sea level contribution from land water storage changes (m).

    # --------------------
    # Model Equations
    # --------------------

    function run_timestep(p, v, d, t)

        if gettime(t) < p.first_projection_year

            # Anthropogenic landwater storage values set to 0 during model calibration period (sea level observations removed land water storage contribution).
            v.lws_sea_level[t] = 0.0

        elseif gettime(t) == p.first_projection_year

            # Set initial condition for first model projection year.
            v.lws_sea_level[t] = p.lws₀

        else
            # Add landwater storage values for projection period.
            v.lws_sea_level[t] = v.lws_sea_level[t-1] + p.lws_distribution[t]
            #v.lws_sea_level[t] = v.lws_sea_level[t-1] + rand(p.lws_distribution)
        end
    end
end
