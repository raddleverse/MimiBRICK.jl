module TestModifications

using MimiBRICK
using Test

# The land water storage is probabilistic, so global sea level change, AIS, and
# land water storage won't be exact every time. Check total GMSL against the sum
# of the components, and length of output match, as those won't be dependent on
# any future changes to model components.
testtol = 1E-6

##==============================================================================
## Checking other RCP scenarios, time periods, model configurations

m = MimiBRICK.get_model(rcp_scenario="RCP26", end_year=2300)
run(m)
@test length(m[:global_sea_level,:sea_level_rise]) == 451
tot = m[:antarctic_icesheet,:ais_sea_level][end] + m[:glaciers_small_icecaps,:gsic_sea_level][end] +
      m[:greenland_icesheet,:greenland_sea_level][end] + m[:thermal_expansion,:te_sea_level][end] +
      m[:landwater_storage,:lws_sea_level][end]
@test m[:global_sea_level,:sea_level_rise][end] ≈ tot atol = testtol

m = MimiBRICK.get_model(rcp_scenario="RCP45", start_year=1900, end_year=2300)
run(m)
@test length(m[:global_sea_level,:sea_level_rise]) == 401
tot = m[:antarctic_icesheet,:ais_sea_level][end] + m[:glaciers_small_icecaps,:gsic_sea_level][end] +
      m[:greenland_icesheet,:greenland_sea_level][end] + m[:thermal_expansion,:te_sea_level][end] +
      m[:landwater_storage,:lws_sea_level][end]
@test m[:global_sea_level,:sea_level_rise][end] ≈ tot atol = testtol

m = MimiBRICK.get_model(rcp_scenario="RCP60", start_year=1950, end_year=2300)
run(m)
@test length(m[:global_sea_level,:sea_level_rise]) == 351
tot = m[:antarctic_icesheet,:ais_sea_level][end] + m[:glaciers_small_icecaps,:gsic_sea_level][end] +
      m[:greenland_icesheet,:greenland_sea_level][end] + m[:thermal_expansion,:te_sea_level][end] +
      m[:landwater_storage,:lws_sea_level][end]
@test m[:global_sea_level,:sea_level_rise][end] ≈ tot atol = testtol

m = MimiBRICK.get_model(rcp_scenario="RCP85", start_year=2000, end_year=2300)
run(m)
@test length(m[:global_sea_level,:sea_level_rise]) == 301
tot = m[:antarctic_icesheet,:ais_sea_level][end] + m[:glaciers_small_icecaps,:gsic_sea_level][end] +
      m[:greenland_icesheet,:greenland_sea_level][end] + m[:thermal_expansion,:te_sea_level][end] +
      m[:landwater_storage,:lws_sea_level][end]
@test m[:global_sea_level,:sea_level_rise][end] ≈ tot atol = testtol

end