# Include all necessary calibration functions

include("calibration_helper_functions.jl")

include("run_historic_models/run_brick_historic_climate.jl")
include("run_historic_models/run_doeclimbrick_historic_climate.jl")
include("run_historic_models/run_sneasybrick_historic_climate.jl")

include("create_log_posteriors/create_log_posterior_brick.jl")
include("create_log_posteriors/create_log_posterior_doeclimbrick.jl")
include("create_log_posteriors/create_log_posterior_sneasybrick.jl")

include("calibration.jl")
