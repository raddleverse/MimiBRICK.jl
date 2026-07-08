using MimiBRICK
using Turing
using DynamicPPL
using ADTypes
using CSVFiles
using DataFrames
using LinearAlgebra

##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
## This script tests whether the SNEASY-BRICK Turing model (`construct_sneasybrick_turing_model`,
## defined in src/calibration/create_log_posteriors/create_turing_model_sneasybrick.jl) can be
## calibrated with Turing's NUTS sampler.
##
## Short answer: not with NUTS's default (ForwardDiff) autodiff backend. See "Result" below.
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------

calibration_start_year = 1850
calibration_end_year   = 2017

run_mymodel! = MimiBRICK.construct_run_sneasybrick(calibration_start_year, calibration_end_year)

model, empirical_R, eta = MimiBRICK.construct_sneasybrick_turing_model(run_mymodel!;
    model_start_year = calibration_start_year, calibration_end_year = calibration_end_year,
    joint_antarctic_prior = false, uniform_ECS = false)

##------------------------------------------------------------------------------
## NUTS (and any other gradient- or Hamiltonian-based sampler) needs a starting point where
## both the log density AND its gradient are finite. Raw draws from the priors essentially
## never give a physically sensible SNEASY-BRICK run (this is a pre-existing property of this
## model, not something introduced by the Turing conversion — `run_calibration`'s own
## `start_from_priors=true` option is explicitly unimplemented for the same reason). So we seed
## the chain from the same known-good starting values the package's RAM-sampler calibration uses.
##------------------------------------------------------------------------------

calibration_data_dir = joinpath(pkgdir(MimiBRICK), "data", "calibration_data")
initvals = DataFrame(load(joinpath(calibration_data_dir, "calibration_initial_values_sneasybrick_04Jun2026.csv"))).parameter_values

# Same order as `p` in `sneasybrick_model` / the initial-values CSV (36 scalar parameters, then
# the 15 Antarctic ice sheet parameters as a block).
scalar_param_names = (
    :σ_temperature, :σ_ocean_heat, :σ_glaciers, :σ_greenland, :σ_antarctic, :σ_gmsl, :σ²_white_noise_CO₂,
    :ρ_temperature, :ρ_ocean_heat, :ρ_glaciers, :ρ_greenland, :ρ_antarctic, :ρ_gmsl, :α₀_CO₂,
    :CO₂_0, :N₂O_0, :temperature_0, :ocean_heat_0, :thermal_s₀, :greenland_v₀, :glaciers_v₀, :glaciers_s₀, :antarctic_s₀,
    :Q10, :CO₂_fertilization, :CO₂_diffusivity,
    :heat_diffusivity, :rf_scale_aerosol, :ECS,
    :thermal_α,
    :greenland_a, :greenland_b, :greenland_α, :greenland_β,
    :glaciers_β₀, :glaciers_n,
)

initial_params = merge(
    NamedTuple(zip(scalar_param_names, initvals[1:36])),
    (antarctic_params = Vector{Float64}(initvals[37:51]), R = empirical_R),
)

##------------------------------------------------------------------------------
## Attempt 1: NUTS with the default (ForwardDiff) autodiff backend.
##------------------------------------------------------------------------------

println("Attempting NUTS with ForwardDiff (Turing's default AD backend)...")
try
    chain = sample(model, NUTS(5, 0.65), 5; initial_params = DynamicPPL.InitFromParams(initial_params, nothing))
    println(chain)
    println("NUTS (ForwardDiff) succeeded.")
catch e
    println()
    println("NUTS (ForwardDiff) FAILED, as expected. Root cause:")
    println()
    println("  SNEASY-BRICK is run through a single `Mimi.ModelInstance` that is built once, up")
    println("  front, in `construct_run_sneasybrick` (run_sneasybrick_historic_climate.jl) and")
    println("  reused (via closure) for every posterior evaluation. Mimi stores each of that")
    println("  instance's scalar parameters in a concretely-typed `Mimi.ScalarModelParameter{Float64}`")
    println("  (the type is fixed once, at `Mimi.build()` time). NUTS's default AD backend")
    println("  (ForwardDiff) needs to push `ForwardDiff.Dual` numbers through the model to get a")
    println("  gradient, but `update_param!` can't store a `Dual` into that `Float64`-typed slot:")
    println("  it fails with `MethodError: no method matching Float64(::ForwardDiff.Dual{...})`")
    println("  inside `Mimi.jl` itself (connections.jl, `setproperty!` on `ScalarModelParameter`).")
    println()
    println("  This is a Mimi.jl architecture limitation (parameter storage type is fixed at build")
    println("  time), not a bug in this repo's SNEASY-BRICK components — those already declare their")
    println("  Mimi parameters generically (`Parameter()`, no explicit `Float64`), and the calibration")
    println("  wrapper code (create_turing_model_sneasybrick.jl, calibration_helper_functions.jl,")
    println("  run_sneasybrick_historic_climate.jl) has been generalized to preserve whatever numeric")
    println("  type flows through `p` rather than truncating it to `Float64`. It's not enough on its")
    println("  own: the fix would have to happen inside Mimi.jl (or by rebuilding a fresh")
    println("  `Mimi.ModelInstance` — a much heavier per-evaluation cost — for every distinct AD tag).")
    println()
    println("Exception detail:")
    showerror(stdout, e)
    println()
end

##------------------------------------------------------------------------------
## Attempt 2 (commented out by default): NUTS with a finite-difference AD backend instead.
## `AutoFiniteDiff()` never constructs a `Dual` — it just repeatedly evaluates the model at
## perturbed `Float64` points — so it sidesteps the `MethodError` above entirely. But each
## gradient then costs roughly (#parameters + 1) full SNEASY-BRICK runs (~75 here), and NUTS
## needs several gradient evaluations per sample, so this is impractically slow for real use.
## Left here for the curious; expect single-digit *minutes* per sample, not seconds.
##------------------------------------------------------------------------------

# println("Attempting NUTS with AutoFiniteDiff (slow, but avoids the ForwardDiff/Mimi conflict)...")
# chain = sample(model, NUTS(2, 0.65; adtype = AutoFiniteDiff()), 3;
#                initial_params = DynamicPPL.InitFromParams(initial_params, nothing))
# println(chain)

##------------------------------------------------------------------------------
## What actually works today: a gradient-free Turing sampler, e.g. MH. See
## src/calibration/create_log_posteriors/create_turing_model_sneasybrick.jl for the model
## definition; `sample(model, MH(), N)` calibrates it without needing any of the above.
##------------------------------------------------------------------------------

println()
println("Falling back to MH (gradient-free, works today) for a short demonstration run...")
chain_mh = sample(model, MH(), 5; initial_params = DynamicPPL.InitFromParams(initial_params, nothing))
println(chain_mh)
