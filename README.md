![](https://github.com/raddleverse/MimiBRICK.jl/workflows/Run%20CI%20on%20master/badge.svg)
[![codecov](https://codecov.io/gh/raddleverse/MimiBRICK.jl/branch/master/graph/badge.svg?token=H7SJB47W5V)](https://codecov.io/gh/raddleverse/MimiBRICK.jl)

# MimiBRICK.jl

This repository will have the latest "standard" BRICK version and codes demonstrating its calibration in three configurations: (i) alone (forced by temperature and ocean heat), (ii) coupled to DOECLIM, and (iii) coupled to SNEASY. Standard calibration output will be provided as well as examples and scripts for downscaling the projections to a 1-degree grid for local sea-level rise.

## Recommended Use

* If you would like to use previously published BRICK results as part of a new analysis, but do not necessarily want to re-run the model, then we recommend to go to the [**Introduction-and-Library** repository](https://github.com/raddleverse/Introduction-and-Library). There, you will find links for each published BRICK study, including links for the calibrated model parameter data sets and the calibrated model projections for sea level, temperature, and any other relevant outputs.
* If you would like to run the model yourself, then you are in the right place!
  * Load the MimiBRICK package (or if you with to edit the package either fork the master branch from this repository or download the zipped file of codes)
  * Do your analysis
  * Please let us know of any new results that should be incorporated into the [Introduction-and-Library repository](https://github.com/raddleverse/Introduction-and-Library)! See the README.md file in that repository for more information.
  * Also please let us know of any model modifications and/or bug fixes that might usefully be incorporated into the main BRICK codes. Creating an "Issue" here is a great way to do that (top horizontal menu bar in the GitHub browser).

### License

Copyright 2022 Tony Wong, (many others - TODO)

This file is part of MimiBRICK.jl (Building blocks for Relevant Ice and Climate Knowledge). MimiBRICK.jl is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

MimiBRICK.jl is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with MimiBRICK.jl (`LICENSE.md)`). If not, see http://www.gnu.org/licenses/.

## How To Install MimiBRICK

This code was created using [Julia v1.6](https://julialang.org/downloads/) and requires several Julia packages. It is recommended that you use Julia v1.6 (or later). Julia may be downloaded from http://julialang.org/downloads/.

(1) Run the following line to install the Mimi implementation of BRICK:

```julia
] add MimiBRICK
```

(2) To exit back to Julia, hit the `backspace` key.

## Running baseline cases with default parameters and unit tests

If you would like to take a look at the unit tests run with Continuous Integration for this package, feel free to take a look at the `runtests.jl` file, and the separate testing scripts it calls [here](https://github.com/raddleverse/MimiBRICK.jl/tree/master/test).

### BRICK standalone (with temperature and ocean heat uptake exogenous forcing)

This is the first test that is done in `test/runtests/test_default.jl`. Since it does not require DOECLIM or SNEASY, you can run BRICK using temperature and ocean heat uptake forcing data by running in the Julia console:
```julia
using MimiBRICK
m = MimiBRICK.get_model()
run(m)
```

You can plot the output fields in the model object `m` using (for example) the `Plots` Julia plotting package. First, let's grab the years over which the model was run. This is a dimension in the model. We can retrieve it by using the `dim_keys()` function, from the `Mimi` package.
```julia
using Mimi
years = dim_keys(m, :time)
```

Then we can load the `Plots` package and make a figure of the global mean sea-level change. Note that the first argument into the `m` object specifies the component of our model, and the second argument specifies the field. Here, we are grabbing the `sea_level_rise` field from the `global_sea_level` component.
```julia
using Plots
plot(years, m[:global_sea_level, :sea_level_rise])
```

Mimi also offers an explorer window to check these model output fields out. To use this, we need to load the `Mimi` package (if you haven't already).
```julia
using Mimi
```

Then, we can open the explorer.
```julia
explore(m)
```

This should open a window labeled "Mimi Explorer Window". On the left, there should be two vertically-stacked boxes. The top box is labeled "Components" and the bottom box is labeled "Data". To view some of the model output from our out-of-box BRICK simulation, you must first pick one of the Components from the top box, and then a Variable field out of the bottom box. For example, to view the global mean sea level model projections, select `global_sea_level` from the Components box, and `global_sea_level : sea_level_rise` from the Variables box. A plot of this variable should appear in the right column of plot boxes. The top box is static; the bottom box you can interact with to zoom in on different portions of the figure.

More information about exploring Mimi model results can be found in the [Mimi Framework How-To guides online](https://www.mimiframework.org/Mimi.jl/stable/howto/howto_2/).

### DOECLIM-BRICK (with radiative forcing)

A simulation using DOECLIM to model temperature and ocean heat uptake, coupled to BRICK for sea-level rise can be constructed and run using the `create_brick_doeclim()` function. This function uses the `MimiSNEASY` package, which includes the DOECLIM model. Here we construct the model object `m` using the `create_brick_doeclim()` constructor, and run the model using the standard Mimi `run(m)`.

```julia
m = MimiBRICK.create_brick_doeclim()
run(m)
```

These are using the default arguments in the model constructor for the DOECLIM-BRICK model:
* `rcp_scenario = "RCP85"` - using Representative Concentration Pathway 8.5 as a default; other options include `RCP26`, `RCP45`, and `RCP60`
* `start_year = 1850` - starting year of the model simulation
* `end_year = 2020` - ending year of the model simulation

So, if you wanted to instead run DOECLIM-BRICK using RCP 6.0 from 1800 to 2100, you could run:
```julia
m = MimiBRICK.create_brick_doeclim(rcp_scenario="RCP60", start_year=1800, end_year=2100)
run(m)
```

And of course, you can use `explore(m)` to check out the model outputs attached to the model object `m` to verify that we have in fact changed RCP scenario and time periods. (Reminder: you'll need to have loaded the `Mimi` package using `using Mimi` to access the `explore()` function.)

### SNEASY-BRICK (with emissions forcing)

Running a coupled model using SNEASY and BRICK proceeds in the same way as DOECLIM-BRICK. We can use the constructor for the coupled SNEASY-BRICK model, create, then run the coupled model.
```julia
m = MimiBRICK.create_sneasy_brick()
run(m)
```

The `create_sneasy_brick()` constructor has the same arguments as the DOECLIM-BRICK constructor, so you can change the RCP scenario and the time period.

## Running the model calibration

**Warning: expert users only!** All others - this is the calibration that is performed using Markov chain Monte Carlo. It leads to the parameter sub-samples that are used for analysis, described below. You do not necessarily need to mess around with this part of the code. If you blindly start running `examples/calibration_driver.jl` out of the box, it will take a **long** time (8-20 hours, probably).

The calibration that is done here follows the same procedure as outlined in [Wong et al. (2017)](https://gmd.copernicus.org/articles/10/2741/2017/) and other works using BRICK. For each of the three main model configurations supported here (BRICK, DOECLIM-BRICK and SNEASY-BRICK), we:
* run a Markov chain Monte Carlo calibration using 20 million iterations
* remove at least 1 million iterations from the beginning of the Markov chain for burn-in
  * the specific length depends on the model configuration; [Gelman and Rubin (1992)](https://projecteuclid.org/journals/statistical-science/volume-7/issue-4/Inference-from-Iterative-Simulation-Using-Multiple-Sequences/10.1214/ss/1177011136.full) potential scale reduction factor is checked < 1.1 for convergence
* subsample 10,000 concomitant parameter sets from the remaining burned-in chain. These samples are used for the hindcast and projections for analysis

This is all done by running the `examples/calibration_driver.jl` script. This script runs the `MimiBRICK.run_calibration` function three times: using `model_config=brick`, `doeclimbrick` and `sneasybrick`. If you want to verify that things are working properly but not wait hours for results, then it is recommended that you try a shorter calibration. This is done in `Calibration_Example.ipynb`, and is achieved by by modifying the arguments for:
* `total_chain_length` - for the three configurations, 1 million iterations typically takes less than an hour. If you are just checking that things are working properly, doing 10,000 would of course be faster, and likely sufficient
* `burnin_length` - this must be less than `total_chain_length`
* `threshold_gr` - if you do a short test calibration, it will yell at you that some of the parameters' potential scale reduction factors are not less than this threshold. You don't need to do anything about it, just letting you know so you don't worry about it.
* `size_subsample` - this must be less than `total_chain_length - burnin_length`

The `run_calibration` function will create a subfolder for the `model_config` in the user-defined `output_dir`. Within that results directory, you will find:
* `parameters_full_chain_(model_config).csv` - the full Markov chain of parameter samples, including the burn-in period
* `mcmc_log_post_(model_config).csv` - the log-posterior scores (numerator from Bayes' theorem) for the full chain of parameter samples
* `parameters_subsample_(model_config).csv` - the parameter values in the sub-sample for analysis
* `log_post_subsample_(model_config).csv` - the log-posterior scores for the sub-sample of parameters for analysis. This is used to determine the maximum _a posteriori_ simulation
* `proposal_covariance_matrix_(model_config).csv` - the final proposal covariance matrix for the adaptive proposals. If you use this and the final sample of parameters from `parameters_full_chain_(model_config).csv`, you can restart the Markov chain calibration. This and the last iteration of the Markov chain are both saved under the `calibration_data/from_calibration_chains` subdirectory.
* `mcmc_acceptance_rate_(model_config).csv` - the acceptance rate from the MCMC algorithm. Should be about 0.23 for the numbers of parameters (dimension) that we're dealing with here.

Note that calibrations of 20 million iterations will take multiple hours to complete.
* For BRICK on its own, this took about 8 hours on a standard desktop workstation (ca. 2020)
* For DOECLIM-BRICK and SNEASY-BRICK, this will take closer to 15 hours or so (using that same machine)

## Running the model hindcasts

This is done for the hindcast period 1850-2017 by using the `MimiBRICK.run_hindcast` function, using `model_config=brick`, `doeclimbrick` and `sneasybrick` (three times). For the hindcast, no RCP scenario needs to be specified, because all of them follow historical radiative forcing/emissions trends up to 2005.

The standard set of parameters that are being used for the hindcast and projection simulations are the sub-sample of 10,000 from the MCMC calibration described above (`parameters_subsample_(model_config).csv`). If you have a different parameter file that you want to run the hindcasts under, you will want to modify the section of `run_hindcast.jl` titled `Set paths for results files` (line 41).

This script will add model configuration-specific directory that was constructed above (or came with the model codes). It will create a sub-directory called `hindcast_csv` which will be populated with CSV files that include the simulated hindcasts of the model output fields. Each of these names is appended with `model_config` (`brick`, `doeclimbrick`, or `sneasybrick`) and contains one hindcast simulation for each of the sets of parameters in the sub-sample for analysis. Rows correspond to different years (1850-2017 be default) and columns each correspond to different ensemble members.
* `hindcast_antarctic_(model_config).csv` - contribution to global mean sea-level change from the Antarctic ice sheet (meters)
* `hindcast_greenland_(model_config).csv` - contribution to global mean sea-level change from the Greenland ice sheet (meters)
* `hindcast_glaciers_(model_config).csv` - contribution to global mean sea-level change from glaciers and ice caps (meters)
* `hindcast_landwater_storage_sl_(model_config).csv` - contribution to sea-level change from land water storage (meters)
* `hindcast_gmsl_(model_config).csv` - total global mean sea-level change (meters)
* `hindcast_ocean_heat_(model_config).csv` - (DOECLIM- or SNEASY-BRICK only)
* `hindcast_temperature_(model_config).csv` - (DOECLIM- or SNEASY-BRICK only)
* `hindcast_oceanco2_(model_config).csv` - (SNEASY-BRICK only)
* `hindcast_co2_(model_config).csv` - (SNEASY-BRICK only)
* `hindcast_MAP_(model_config).csv` - all of the hindcast time series for the maximum _a posteriori_ set of parameters

## Running the model projections under different RCP scenarios

This is done for the period 1850-2300 (but can be modified to any period between 1765 and 2300) by using the `MimiBRICK.run_projections` function, using `model_config=brick`, `doeclimbrick` or `sneasybrick` and `rcp_scenario="RCP26"`, `"RCP45"`, `"RCP60"`, or `"RCP85"`. Note that the RCP scenario forcing files are all the same until 2005, and the provided stand-alone BRICK temperature and ocean heat forcing files cover the period 1850-2300.

This script will add to the model configuration-specific directory that was constructed above (or came with the model codes). It will create a sub-directory called `projections_csv`, and a sub-directory within there that is specific to each RCP scenario used will be created. The projections files are analogous to the hindcast files that are generated, and will populate the `projections_csv/[RCP scenario]` directory.

## Creating forcing files for stand-alone BRICK

The forcing files for DOECLIM-BRICK (radiative forcing) and SNEASY-BRICK (emissions) are taken from the RCP database here (https://tntcat.iiasa.ac.at/RcpDb/dsd?Action=htmlpage&page=download) and the data repository of Malte Meinshausen here (http://www.pik-potsdam.de/~mmalte/rcps/data/).

For stand-alone BRICK, which requires temperature and ocean heat time series as forcing data, we use the time series for temperature and ocean heat uptake from the maximum _a posteriori_ simulations from the SNEASY-BRICK simulation ensembles described above. This is done in the `calibration/sneasy_make_brick_forcing.jl` script. This script creates the following files, where `xx` denotes the RCP scenario (`26`, `45`, `60`, or `85`), `yyyy` denotes the starting year of the forcing, and `YYYY` denotes the ending year of the forcing.
* `data/model_data/sneasy_oceanheat_RCPxx_yyyy_YYYY.csv`
* `data/model_data/sneasy_temperature_RCPxx_yyyy_YYYY.csv`

## Generating projections of local mean sea-level change

The `MimiBRICK.downscale_brick` function downscales the BRICK global sea level projections to local. This uses the sea-level "fingerprints" of [Slangen et al. (2014)](https://link.springer.com/article/10.1007/s10584-014-1080-9). The downscaling routine will automatically create a subdirectory in the output directory called `localslr`. In this subdirectory, the routine will save an output file with the downscaled local mean sea level change model output.

This routine will downscale either a full ensemble of BRICK model simulations or just the maximum a posteriori model simulation to a specific latitude and longitude point. These are provided by the user as `lat` (degrees north) and `lon` (degrees east). Other needed function arguments include:
* `results_dir` - (String) the directory holding model outputs
* `model_config` - (String) one of `"brick"`, `"doeclimbrick"`, or `"sneasybrick"`. Only the BRICK projections are being downscaled (no CO2 or temperature, for example), but the `downscale_brick` function will find the relevant input data and tag the output files appropriately based on the `model_config` setting.
* `proj_or_hind` - (String) one of `"proj"` (projections) or `"hind"` (hindcast). They're treated similarly when running the model, but this helps for finding the output files in the `results_dir` directories.
* `rcp_scenario` - (String) one of `"RCP26"`, `"RCP45"`, `"RCP60"`, or `"RCP85"`. If running a hindcast, this does not matter.
* `ensemble_or_map` - (String) one of `"ensemble"` or `"map"`. If `"ensemble"`, then will downscale the full BRICK ensemble that matches the provided `model_config`, `proj_or_hind`, and `rcp_scenario` settings. If `"map"`, will only downscale the maximum a posteriori simulation.

An example for use can be found in `test/test_downscaling.jl` as well as `examples/Calibration_Example.jl`. This example includes a few specifications, for example in one the BRICK (standalone model) maximum a posteriori sea-level rise projection under RCP8.5 is downscaled for New York City using the following settings:

```julia
lat = 40.7128 # deg N
lon = 360-74.0060 # 74.0060 deg W
model_config = "brick"
proj_or_hind = "proj"
rcp_scenario = "RCP85"
ensemble_or_map = "map"
```

The following line of code performs the actual downscaling and saves the output files to the appropriate `results_dir` directory.
```julia
years, lsl = downscale_brick(lon=lon, lat=lat, results_dir=results_dir, proj_or_hind=proj_or_hind, ensemble_or_map=ensemble_or_map, model_config=model_config, rcp_scenario=rcp_scenario)
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
