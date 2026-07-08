using MimiBRICK

##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
## Runs a per-block Robust Adaptive Metropolis (RAM; Vihola 2012) calibration of the SNEASY-BRICK
## Turing model (`MimiBRICK.sneasybrick_model`, see `create_turing_model_sneasybrick.jl`). See
## `src/calibration/turing_ram_calibration.jl` for the sampler itself and the design rationale
## (why a per-block *adaptive* proposal, rather than a fixed-covariance Gibbs/MH block, is needed
## here).
##
## Override the sweep count either on the command line
## (`julia --project=. examples/turing_mh_ram_calibration.jl 500`) or by calling
## `MimiBRICK.run_per_block_ram_calibration(N; output_dir=...)` directly with any other keyword
## arguments (e.g. `reference_posterior_path` to warm-start from a previous calibration's posterior
## subsample instead of the bundled default).
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------

output_dir = joinpath(@__DIR__, "results")

n_samples = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2_000

df, elapsed, acceptance_rates = MimiBRICK.run_per_block_ram_calibration(n_samples; output_dir = output_dir)
