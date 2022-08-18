# MimiBRICK.jl

```@contents
```

### Functions for creating typical coupled model instances beyond stand-alone BRICK.

```@docs

MimiBRICK.create_brick_doeclim(;rcp_scenario::String = "RCP85", start_year::Int=1850, end_year::Int=2020)

MimiBRICK.create_sneasy_brick(;rcp_scenario::String = "RCP85", start_year::Int=1850, end_year::Int=2020)
```

### Functions for creating prior and posterior distributions for BRICK and typical coupled model configuration (including likelihoods)

```@docs
MimiBRICK.construct_brick_log_prior(joint_antarctic_prior::Bool; calibration_data_dir::Union{String, Nothing} = nothing)

MimiBRICK.construct_brick_log_posterior(f_run_model!; model_start_year::Int=1852, calibration_end_year::Int=2017, joint_antarctic_prior::Bool=false)

MimiBRICK.construct_doeclimbrick_log_prior(joint_antarctic_prior::Bool, uniform_ECS::Bool; calibration_data_dir::Union{String, Nothing} = nothing)

MimiBRICK.construct_doeclimbrick_log_posterior(f_run_model!; model_start_year::Int=1850, calibration_end_year::Int=2017, joint_antarctic_prior::Bool=false, uniform_ECS::Bool=false)

MimiBRICK.construct_sneasybrick_log_prior(joint_antarctic_prior::Bool, uniform_ECS::Bool; calibration_data_dir::Union{String, Nothing} = nothing)

MimiBRICK.construct_sneasybrick_log_posterior(f_run_model!; model_start_year::Int=1850, calibration_end_year::Int=2017, joint_antarctic_prior::Bool=false, uniform_ECS::Bool=false)
```

### Calibration functionality

```@docs
MimiBRICK.load_calibration_data(model_start_year::Int, last_calibration_year::Int; last_sea_level_norm_year::Int=1990, calibration_data_dir::Union{Nothing, String} = nothing)

MimiBRICK.hetero_logl_ar1(residuals::Array{Float64,1}, σ::Float64, ρ::Float64, ϵ::Array{Union{Float64, Missings.Missing},1})

MimiBRICK.hetero_logl_car1(residuals::Array{Float64,1}, indices::Array{Int64,1}, σ²_white_noise::Float64, α₀::Float64, ϵ::Array{Union{Float64, Missings.Missing},1})

MimiBRICK.calculate_trends(model_output::Array{Float64,1}, obs_trends::DataFrame, start_year::Int, end_year::Int)

MimiBRICK.truncated_kernel(data, lower_bound, upper_bound)

MimiBRICK.simulate_ar1_noise(n::Int, σ::Float64, ρ::Float64, ϵ::Array{Float64,1})

MimiBRICK.simulate_car1_noise(n, α₀, σ²_white_noise, ϵ)

MimiBRICK.replicate_errors(start_year::Int, end_year::Int, error_data)

MimiBRICK.run_calibration(;   output_dir::String,
                      model_config="brick",
                      calibration_start_year=1850,
                      calibration_end_year=2005,
                      total_chain_length=1000,
                      burnin_length=0,
                      threshold_gr=1.1,
                      num_walkers=2,
                      size_subsample=1000,
                      start_from_priors=false,
                      calibration_data_dir::Union{String, Nothing} = nothing
                  )

MimiBRICK.run_hindcast(; output_dir::String,
                  model_config::String = "brick",
                  start_year::Int = 1850,
                  end_year = 2017,
              )

MimiBRICK.run_projections(; output_dir::String,
                      model_config::String = "brick",
                      rcp_scenario::String = "RCP85",
                      start_year::Int = 1850,
                      end_year = 2300,
                  )
```
