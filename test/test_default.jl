module TestDefault

using MimiBRICK
using Test

# The land water storage is probabilistic, so global sea level change, AIS, and
# land water storage won't be exact every time. Check total GMSL against the sum
# of the components, and length of output match, as those won't be dependent on
# any future changes to model components.
testtol = 1E-6

##==============================================================================
## With default parameters

## BRICK standalone
m = MimiBRICK.get_model()
run(m)
## Run tests for the out-of-box model with default parameters and forcing
@test length(m[:global_sea_level,:sea_level_rise]) == 171
tot = m[:antarctic_icesheet,:ais_sea_level][end] + m[:glaciers_small_icecaps,:gsic_sea_level][end] +
      m[:greenland_icesheet,:greenland_sea_level][end] + m[:thermal_expansion,:te_sea_level][end] +
      m[:landwater_storage,:lws_sea_level][end]
@test m[:global_sea_level,:sea_level_rise][end] ≈ tot atol = 0.000001

## DOECLIM+BRICK
m = MimiBRICK.create_brick_doeclim()
run(m)
@test length(m[:global_sea_level,:sea_level_rise]) == 171
tot = m[:antarctic_icesheet,:ais_sea_level][end] + m[:glaciers_small_icecaps,:gsic_sea_level][end] +
      m[:greenland_icesheet,:greenland_sea_level][end] + m[:thermal_expansion,:te_sea_level][end] +
      m[:landwater_storage,:lws_sea_level][end]
@test m[:global_sea_level,:sea_level_rise][end] ≈ tot atol = testtol

## SNEASY+BRICK
m = MimiBRICK.create_sneasy_brick()
run(m)
println(m[:ccm,:atmco2][end])
@test length(m[:global_sea_level,:sea_level_rise]) == 171
tot = m[:antarctic_icesheet,:ais_sea_level][end] + m[:glaciers_small_icecaps,:gsic_sea_level][end] +
      m[:greenland_icesheet,:greenland_sea_level][end] + m[:thermal_expansion,:te_sea_level][end] +
      m[:landwater_storage,:lws_sea_level][end]
@test m[:global_sea_level,:sea_level_rise][end] ≈ tot atol = testtol

end