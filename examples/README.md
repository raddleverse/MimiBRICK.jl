# Calibration, Hindcast/Projection Model Simulations, and Downscaling

## Running the model calibration

**Warning: expert users only!** All others - this is the calibration that is performed using Markov chain Monte Carlo. It leads to the parameter sub-samples that are used for analysis, described below. You do not necessarily need to mess around with this part of the code. If you blindly start running `examples/calibration_driver.jl` out of the box, it will take a long time (about 10 hours, probably, depending on the computer you're running it on).

The calibration that is done here follows the same procedure as outlined in [Wong et al. (2017)](https://gmd.copernicus.org/articles/10/2741/2017/) and other works using BRICK. For each of the three main model configurations supported here (BRICK, DOECLIM-BRICK and SNEASY-BRICK), we:
* run an initial Markov chain Monte Carlo sampler for 15 million iterations, to improve initial conditions for the "production" sampling later
* run a Markov chain Monte Carlo sampler using 15 million iterations
* remove 5 million iterations from the beginning of the Markov chain for burn-in
  * [Gelman and Rubin (1992)](https://projecteuclid.org/journals/statistical-science/volume-7/issue-4/Inference-from-Iterative-Simulation-Using-Multiple-Sequences/10.1214/ss/1177011136.full) potential scale reduction factor is checked < 1.1 for convergence
* subsample 10,000 concomitant parameter sets from the remaining burned-in chain. These samples are used for the hindcast and projections for analysis

This is all done by running the `examples/calibration_driver.jl` script. This script runs the `MimiBRICK.run_calibration` function three times: using `model_config=brick`, `doeclimbrick` and `sneasybrick`. If you want to verify that things are working properly but not wait hours for results, then it is recommended that you try a shorter calibration. This is done in `Calibration_Example.ipynb`, and is achieved by by modifying the arguments for:
* `total_chain_length` - for the three configurations, 1 million iterations typically takes less than an hour. If you are just checking that things are working properly, doing 10,000 would of course be faster, and likely sufficient
* `burnin_length` - this must be less than `total_chain_length`
* `threshold_gr` - if you do a short test calibration, it will yell at you that some of the parameters' potential scale reduction factors are not less than this threshold. You don't need to do anything about it, just letting you know so you don't worry about it.
* `size_subsample` - this must be less than `total_chain_length - burnin_length`
* `gmsl_data` - `:wa` (Wang et al. 2024, the default) or `:cw` (Church & White 2011). This option applies to all three model configurations.
* `calibration_targets` - `:standard` for MimiBRICK's usual mixture of calibration datasets, or `:mengel_ext` for Marcus Sarofim's extended Mengel target table. The `:mengel_ext` option uses Frederikse-based AIS, glacier, Greenland, steric, and LWS series with Dangendorf/NOAA total GMSL, all referenced to 1995–2005. It requires `glacier_model=:mengel`; use `calibration_end_year=2026` to include every available target year.

When `calibration_targets=:mengel_ext`, calibration and projection files use
the suffix `_data-mengel-ext` instead of `_gmsl-wa` or `_gmsl-cw`.

The `run_calibration` function will create a subfolder for the `model_config` in the user-defined `output_dir`. Within that results directory, you will find:
* `parameters_full_chain_(model_config)_gmsl-(gmsl_data).csv` - the full Markov chain of parameter samples, including the burn-in period
* `mcmc_log_post_(model_config)_gmsl-(gmsl_data).csv` - the log-posterior scores (numerator from Bayes' theorem) for the full chain of parameter samples
* `parameters_subsample_(model_config)_gmsl-(gmsl_data).csv` - the parameter values in the sub-sample for analysis
* `log_post_subsample_(model_config)_gmsl-(gmsl_data).csv` - the log-posterior scores for the sub-sample of parameters for analysis. This is used to determine the maximum _a posteriori_ simulation
* `proposal_covariance_matrix_(model_config)_gmsl-(gmsl_data).csv` - the final proposal covariance matrix for the adaptive proposals. If you use this and the final sample of parameters from `parameters_full_chain_(model_config)_gmsl-(gmsl_data).csv`, you can restart the Markov chain calibration. This and the last iteration of the Markov chain are both saved under the `calibration_data/from_calibration_chains` subdirectory.
* `mcmc_acceptance_rate_(model_config)_gmsl-(gmsl_data).csv` - the acceptance rate from the MCMC algorithm. Should be about 0.23 for the numbers of parameters (dimension) that we're dealing with here.

Note that calibrations of 15 million iterations will take multiple hours to complete.
* For BRICK on its own, this took about 2-3 hours on a standard desktop workstation (ca. 2022)
* For DOECLIM-BRICK and SNEASY-BRICK, this will take closer to 4 hours or so (using that same machine)

## Running the model projections under different SSP-RCP scenarios

This is done for the period 1850-2300 (but can be modified to any period between 1750 and 2300) by using the `MimiBRICK.run_projections` function, using `model_config=brick`, `doeclimbrick` or `sneasybrick` and `ssprcp_scenario` set to one of `ssp119, ssp126, ssp245, ssp370, ssp460, ssp585, ssp534-over`. Note that the SSP-RCP scenario forcing files are all the same until 2015, and the provided stand-alone BRICK temperature and ocean heat forcing files cover the period 1850-2300. These are the MAP simulation results from SNEASY-BRICK over that time period.

This script will add to the model configuration-specific directory that was constructed above (or came with the model codes). It will create a sub-directory called `projections_csv`, and a sub-directory within there that is specific to each SSP-RCP scenario used will be created. The projections files are analogous to the hindcast files that are generated, and will populate the `projections_csv/[scenario]` directory.

## Creating forcing files for stand-alone BRICK

The forcing files for DOECLIM-BRICK (radiative forcing) and SNEASY-BRICK (emissions) are taken from the RCMIP database of Nicholls and Lewis (2021; https://doi.org/10.5281/zenodo.4589756).

For stand-alone BRICK, which requires temperature and ocean heat time series as forcing data, we use the time series for temperature and ocean heat uptake from the maximum _a posteriori_ simulations from the SNEASY-BRICK simulation ensembles described above. This is done in the `calibration/sneasy_make_brick_forcing.jl` script. This script creates the following files, where `xx` denotes the scenario, `yyyy` denotes the starting year of the forcing, and `YYYY` denotes the ending year of the forcing.
* `data/model_data/sneasy_oceanheat_sspxx_yyyy_YYYY.csv`
* `data/model_data/sneasy_temperature_sspxx_yyyy_YYYY.csv`

## Generating projections of local mean sea-level change

The `MimiBRICK.downscale_brick` function downscales the BRICK global sea level projections to local. This uses the sea-level "fingerprints" of [Slangen et al. (2014)](https://link.springer.com/article/10.1007/s10584-014-1080-9). The downscaling routine will automatically create a subdirectory in the output directory's `hindcast_csv` or `projections_csv/RCPXX/` subdirectory (depending on specification of `proj_or_hind` argument) called `localslr`. In this subdirectory, the routine will save an output file with the downscaled local mean sea level change model output.

This routine will downscale either a full ensemble of BRICK model simulations or just the maximum a posteriori model simulation to a specific latitude and longitude point. These are provided by the user as `lat` (degrees north) and `lon` (degrees east). Other needed function arguments include:
* `results_dir` - (String) the directory holding model outputs
* `model_config` - (String) one of `"brick"`, `"doeclimbrick"`, or `"sneasybrick"`. Only the BRICK projections are being downscaled (no CO2 or temperature, for example), but the `downscale_brick` function will find the relevant input data and tag the output files appropriately based on the `model_config` setting.
* `proj_or_hind` - (String) one of `"proj"` (projections) or `"hind"` (hindcast). They're treated similarly when running the model, but this helps for finding the output files in the `results_dir` directories.
* `ssprcp_scenario` - (String) a valid SSP-RCP scenario. If running a hindcast, this does not matter.
* `ensemble_or_map` - (String) one of `"ensemble"` or `"map"`. If `"ensemble"`, then will downscale the full BRICK ensemble that matches the provided `model_config`, `proj_or_hind`, and `ssprcp_scenario` settings. If `"map"`, will only downscale the maximum a posteriori simulation.

An example for use can be found in `test/test_downscaling.jl` as well as `examples/Calibration_Example.jl`. This example includes a few specifications, for example in one the BRICK (standalone model) maximum a posteriori sea-level rise projection under RCP8.5 is downscaled for New York City using the following settings:

```julia
lat = 40.7128 # deg N
lon = 360-74.0060 # 74.0060 deg W
model_config = "brick"
proj_or_hind = "proj"
ssprcp_scenario = "ssp245"
ensemble_or_map = "map"
```

The following line of code performs the actual downscaling and saves the output files to the appropriate `results_dir` directory.
```julia
years, lsl = downscale_brick(lon=lon, lat=lat, results_dir=results_dir, proj_or_hind=proj_or_hind, ensemble_or_map=ensemble_or_map, model_config=model_config, ssprcp_scenario=ssprcp_scenario)
```

## License

MimiBRICK is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
