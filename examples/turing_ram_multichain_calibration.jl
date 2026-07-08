using Distributed
using CSV
using CSVFiles
using DataFrames
using MCMCChains
using NetCDF
using Random
using Statistics

##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
## Runs `n_chains` independent per-block RAM calibrations (see
## `src/calibration/turing_ram_calibration.jl`) of the SNEASY-BRICK Turing model in parallel worker
## processes. Each chain builds its own Mimi model instance internally (inside
## `run_per_block_ram_calibration`), so there's no shared mutable state between chains -- separate
## `Distributed` worker processes make that isolation automatic (and side-step any question of
## whether Mimi/Turing internals are thread-safe, which multiple `Threads.@spawn`'d chains sharing
## one Julia process would raise).
##
## After sampling, this:
##   1. Reports Gelman-Rubin Rhat and effective sample size (`MCMCChains.ess_rhat`) across the 4
##      chains for every sampled quantity, to a CSV.
##   2. Saves the pooled post-burn-in posterior (draw x chain, for each of the 51 SNEASY-BRICK
##      parameters plus the log-joint density) to NetCDF.
##   3. Draws a 10,000-member subsample from the pooled posterior and feeds it through MimiBRICK's
##      existing `run_projections` to produce a hindcast and a future projection, consolidating
##      both into NetCDF files (`run_projections` itself writes one CSV per output field; this
##      script reads those back in and writes a single NetCDF per hindcast/projection run).
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------

n_chains             = 4
n_samples            = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 10_000
burnin_fraction       = 0.3    # RAM is still adapting its per-block proposal covariance early on;
                                # matches roughly the 1/3 burn-in fraction `run_calibration` uses
                                # for the non-Turing RAM sampler.
size_subsample        = 10_000
ssprcp_scenario        = "ssp245"
hindcast_start_year   = 1850
hindcast_end_year     = 2017
projection_end_year   = 2300
seed                  = 2024

output_dir = joinpath(@__DIR__, "results", "turing_ram_multichain")
mkpath(output_dir)

##------------------------------------------------------------------------------
## 1. Run `n_chains` chains in parallel worker processes.
##------------------------------------------------------------------------------
addprocs(max(0, n_chains - (nprocs() - 1)); exeflags = "--project=$(Base.active_project())")
@everywhere using MimiBRICK

println("Running $n_chains chains of $n_samples sweeps each on $(nprocs() - 1) worker process(es)...")
chain_results = try
    pmap(1:n_chains) do i
        MimiBRICK.run_per_block_ram_calibration(n_samples;
            output_dir = joinpath(output_dir, "chain_$i"),
            seed = seed + i)
    end
finally
    rmprocs(workers())
end

param_names = names(chain_results[1][1])
n_burnin = round(Int, burnin_fraction * n_samples)
n_kept = n_samples - n_burnin

println("Acceptance rates by chain:")
for (i, (_, elapsed, rates)) in enumerate(chain_results)
    println("  chain $i ($(round(elapsed, digits=1))s): ", rates)
end

##------------------------------------------------------------------------------
## 2. Gelman-Rubin Rhat and effective sample size across chains (post-burn-in), for every sampled
## quantity (51 physical/statistical parameters, 49 R_i_j entries, and the log-joint density).
##------------------------------------------------------------------------------
chain_array = Array{Float64}(undef, n_kept, length(param_names), n_chains)
for (c, (df, _, _)) in enumerate(chain_results)
    chain_array[:, :, c] = Matrix(df[(n_burnin + 1):end, :])
end
chains = MCMCChains.Chains(chain_array, param_names)
ess_df = DataFrame(MCMCChains.ess_rhat(chains))
ess_path = joinpath(output_dir, "ess_rhat_sneasybrick.csv")
CSV.write(ess_path, ess_df)
println("Rhat/ESS across $n_chains chains ($(nrow(ess_df)) quantities): median ESS = ",
        round(median(ess_df.ess), digits = 1), ", max Rhat = ", round(maximum(ess_df.rhat), digits = 4))
println("Saved Rhat/ESS table to $ess_path")

##------------------------------------------------------------------------------
## 3. Save the pooled post-burn-in posterior to NetCDF: one variable per SNEASY-BRICK parameter
## (named via the package's established external naming convention) plus `log_post`, each with
## (draw, chain) dimensions.
##------------------------------------------------------------------------------
posterior_path = joinpath(output_dir, "posteriors_sneasybrick.nc")
isfile(posterior_path) && rm(posterior_path)

internal_names = vcat([string(n) for n in MimiBRICK.SNEASYBRICK_SCALAR_PARAM_NAMES], ["antarctic_params[$i]" for i in 1:15])
external_names = MimiBRICK.SNEASYBRICK_EXTERNAL_PARAM_NAMES

for (internal, external) in zip(internal_names, external_names)
    data = Array{Float64}(undef, n_kept, n_chains)
    for (c, (df, _, _)) in enumerate(chain_results)
        data[:, c] = df[(n_burnin + 1):end, internal]
    end
    nccreate(posterior_path, external, "draw", collect(1:n_kept), "chain", collect(1:n_chains))
    ncwrite(data, posterior_path, external)
end

logpost = Array{Float64}(undef, n_kept, n_chains)
for (c, (df, _, _)) in enumerate(chain_results)
    logpost[:, c] = df[(n_burnin + 1):end, "logjoint"]
end
nccreate(posterior_path, "log_post", "draw", collect(1:n_kept), "chain", collect(1:n_chains))
ncwrite(logpost, posterior_path, "log_post")
ncclose(posterior_path)
println("Saved pooled posterior ($n_kept draws x $n_chains chains) to $posterior_path")

##------------------------------------------------------------------------------
## 4. Pool the post-burn-in draws across chains and draw a `size_subsample`-member subsample (same
## convention `run_calibration` uses) to drive the model-running (hindcast/projection) step -- one
## Mimi run per subsample member, not per raw MCMC draw.
##------------------------------------------------------------------------------
pooled = vcat([df[(n_burnin + 1):end, :] for (df, _, _) in chain_results]...)
n_pooled = nrow(pooled)
size_subsample = min(size_subsample, n_pooled)
Random.seed!(seed)
subsample = pooled[randperm(n_pooled)[1:size_subsample], :]

parameters_out = DataFrame()
for (internal, external) in zip(internal_names, external_names)
    parameters_out[!, external] = subsample[!, string(internal)]
end

println("Drew a $(size_subsample)-member posterior subsample (pooled from $n_pooled retained draws) for hindcasts/projections.")

##------------------------------------------------------------------------------
## 5. Hindcast + future projection via `MimiBRICK.run_projections`, consolidated to NetCDF.
## `run_projections` writes its own CSV per output field into
## `<output_dir>/projections_csv/sneasybrick/<scenario>/`, keyed only by `model_config` and
## `ssprcp_scenario` -- not by start/end year -- so the hindcast and projection each need their own
## `output_dir` to avoid one overwriting the other's CSVs.
##------------------------------------------------------------------------------
const PROJECTION_FIELDS = ("gmsl", "landwater_storage_sl", "glaciers", "greenland", "antarctic",
                           "thermal", "temperature", "ocean_heat", "co2", "oceanco2")

function consolidate_to_netcdf(run_output_dir::String, model_config::String, ssprcp_scenario::String,
                                start_year::Int, end_year::Int, nc_path::String)
    isfile(nc_path) && rm(nc_path)
    years = collect(start_year:end_year)
    field_dir = joinpath(run_output_dir, "projections_csv", model_config, ssprcp_scenario)
    for field in PROJECTION_FIELDS
        csv_path = joinpath(field_dir, "projections_$(field)_$(ssprcp_scenario)_$(model_config).csv")
        isfile(csv_path) || continue
        data = coalesce.(Matrix(DataFrame(load(csv_path))), NaN)
        nccreate(nc_path, field, "year", years, "ensemble_member", collect(1:size(data, 2)))
        ncwrite(data, nc_path, field)
    end
    ncclose(nc_path)
    return nc_path
end

for (label, run_dir, start_year, end_year, nc_name) in (
    ("hindcast",   joinpath(output_dir, "for_hindcast"),   hindcast_start_year, hindcast_end_year,   "hindcast_sneasybrick.nc"),
    ("projection", joinpath(output_dir, "for_projection"), hindcast_start_year, projection_end_year, "projections_sneasybrick_$(ssprcp_scenario).nc"),
)
    mkpath(run_dir)
    save(joinpath(run_dir, "parameters_subsample_sneasybrick.csv"), parameters_out)
    save(joinpath(run_dir, "log_post_subsample_sneasybrick.csv"), DataFrame(log_post = subsample.logjoint))

    println("Running $label ($start_year-$end_year, $ssprcp_scenario)...")
    MimiBRICK.run_projections(output_dir = run_dir, model_config = "sneasybrick", ssprcp_scenario = ssprcp_scenario,
        start_year = start_year, end_year = end_year)

    nc_path = consolidate_to_netcdf(run_dir, "sneasybrick", ssprcp_scenario, start_year, end_year, joinpath(output_dir, nc_name))
    println("Saved $label to $nc_path")
end
