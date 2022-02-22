# MimiBRICK.jl

This repository will have the latest "standard" BRICK version and codes demonstrating its calibration in three configurations: (i) alone (forced by temperature and ocean heat), (ii) coupled to DOECLIM, and (iii) coupled to SNEASY. Standard calibration output will be provided as well as examples and scripts for downscaling the projections to a 1-degree grid for local sea-level rise.

## Recommended Use

* If you would like to use previously published BRICK results as part of a new analysis, but do not necessarily want to re-run the model, then we recommend to go to the [**Introduction-and-Library** repository](https://github.com/raddleverse/Introduction-and-Library). There, you will find links for each published BRICK study, including links for the calibrated model parameter data sets and the calibrated model projections for sea level, temperature, and any other relevant outputs.
* If you would like to run the model yourself, then you are in the right place!
  * Either fork the master branch from this repository or download the zipped file of codes.
  * Do your analysis
  * Please let us know of any new results that should be incorporated into the [Introduction-and-Library repository](https://github.com/raddleverse/Introduction-and-Library)! See the README.md file in that repository for more information.
  * Also please let us know of any model modifications and/or bug fixes that might usefully be incorporated into the main BRICK codes. Creating an "Issue" here is a great way to do that (top horizontal menu bar in the GitHub browser).

### License

Copyright 2022 Tony Wong, (many others - TODO)

This file is part of MimiBRICK.jl (Building blocks for Relevant Ice and Climate Knowledge). MimiBRICK.jl is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

MimiBRICK.jl is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with MimiBRICK.jl (`LICENSE.md)`). If not, see http://www.gnu.org/licenses/.


## How To Install Required Packages

This code was created using [Julia v1.6](https://julialang.org/downloads/) and requires several Julia packages.

(1) To install these packages, first enter the package manager by hitting the `]` key in the Julia console. Once in the package manager, run the following code:

```julia
add CSV
add CSVFiles
add DataFrames
add Distributions
add Interpolations
add KernelDensity
add LinearAlgebra
add MCMCDiagnostics
add Mimi
add MimiSNEASY
add NetCDF
add Plots
add RobustAdaptiveMetropolisSampler
add SpecialFunctions
add Statistics
add StatsBase
```

(2) While still in the package manager, run the following line to install the Mimi implementation of SNEASY:

```julia
add https://github.com/anthofflab/MimiSNEASY.jl.git
```

(3) Run the following line to install the Mimi implementation of BRICK:

```julia
add https://github.com/raddleverse/MimiBRICK.jl.git
```

(4) To exit back to Julia, hit the `backspace` key.

## Running baseline cases with default parameters

You'll first want to navigate in your Julia terminal to the `test` directory within this repository. The first three commands in the `runtests.jl` script in that directory will activate the Julia project environment for the `MimiBRICK.jl` codes.

```julia
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()
```

Note that the `".."` in the second command assumes that you are in one of the sub-directories from the main `MimiBRICK.jl` directory. Depending on whether you are working in your own directory system on your own projects, you may decide to modify this.

### BRICK standalone (with temperature and ocean heat uptake exogenous forcing)

This is the first test that is done in `test/runtests.jl`. Since it does not require DOECLIM or SNEASY, you can run BRICK using temperature and ocean heat uptake forcing data by running in the Julia console:
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

A simulation using DOECLIM to model temperature and ocean heat uptake, coupled to BRICK for sea-level rise can be constructed and run. First, we need to load the needed `MimiSNEASY` package, which includes the DOECLIM model. Then we load the `MimiBRICK_DOECLIM.jl` script, which includes a coupled model constructor for DOECLIM-BRICK. Finally, we construct the model object `m` using the `create_brick_doeclim()` constructor, and run the model using the standard Mimi `run(m)`.
```julia
using MimiSNEASY
srcdir = joinpath(@__DIR__, "..", "src")
include(joinpath(srcdir,"MimiBRICK_DOECLIM.jl"))
m = MimiBRICK_DOECLIM.create_brick_doeclim()
run(m)
```

These are using the default arguments in the model constructor for the DOECLIM-BRICK model:
* `rcp_scenario = "RCP85"` - using Representative Concentration Pathway 8.5 as a default; other options include `RCP26`, `RCP45`, and `RCP60`
* `start_year = 1850` - starting year of the model simulation
* `end_year = 2300` - ending year of the model simulation

So, if you wanted to instead run DOECLIM-BRICK using RCP 6.0 from 1800 to 2100, you could run:
```julia
m = MimiBRICK_DOECLIM.create_brick_doeclim(rcp_scenario="RCP60", start_year=1800, end_year=2100)
run(m)
```

And of course, you can use `explore(m)` to check out the model outputs attached to the model object `m` to verify that we have in fact changed RCP scenario and time periods. (Reminder: you'll need to have loaded the `Mimi` package using `using Mimi` to access the `explore()` function.)

### SNEASY-BRICK (with emissions forcing)

Running a coupled model using SNEASY and BRICK proceeds in the same way as DOECLIM-BRICK. Assuming we have already loaded the `MimiSNEASY` package, we can load the constructor for the coupled SNEASY-BRICK model, create, then run the coupled model.
```julia
include(joinpath(srcdir,"create_models","SNEASY_BRICK.jl"))
m = create_sneasy_brick()
run(m)
```

The `create_sneasy_brick()` constructor has the same arguments as the DOECLIM-BRICK constructor, so you can change the RCP scenario and the time period.

## Running the model calibration

**Warning: expert users only!** All others - this is the calibration that is performed using Markov chain Monte Carlo. It leads to the parameter sub-samples that are used for analysis, described below. You do not necessarily need to mess around with this part of the code.

The calibration that is done here follows the same procedure as outlined in [Wong et al. (2017)](https://gmd.copernicus.org/articles/10/2741/2017/) and other works using BRICK. For each of the three main model configurations supported here (BRICK, DOECLIM-BRICK and SNEASY-BRICK), we:
* run a Markov chain Monte Carlo calibration using 20 million iterations
* remove at least 1 million iterations from the beginning of the Markov chain for burn-in
  * the specific length depends on the model configuration; [Gelman and Rubin (1992)](https://projecteuclid.org/journals/statistical-science/volume-7/issue-4/Inference-from-Iterative-Simulation-Using-Multiple-Sequences/10.1214/ss/1177011136.full) potential scale reduction factor is checked < 1.1 for convergence
* subsample 10,000 concomitant parameter sets from the remaining burned-in chain. These samples are used for the hindcast and projections for analysis

This is all done by running the `calibration/calibration.jl` script using `model_config=brick`, `doeclimbrick` and `sneasybrick` (three times).

The script will create a date-stamped directory in `results` specific to this calibration, including the `model_config` and number of Markov chain iterations used. Within that results directory, you will find:
* `parameters_full_chain.csv` - the full Markov chain of parameter samples, including the burn-in period
* `mcmc_log_post.csv` - the log-posterior scores (numerator from Bayes' theorem) for the full chain of parameter samples
* `parameters_subsample.csv` - the parameter values in the sub-sample for analysis
* `log_post_subsample.csv` - the log-posterior scores for the sub-sample of parameters for analysis. This is used to determine the maximum _a posteriori_ simulation
* `proposal_covariance_matrix.csv` - the final proposal covariance matrix for the adaptive proposals. If you use this and the final sample of parameters from `parameters_full_chain.csv`, you can restart the Markov chain calibration. This and the last iteration of the Markov chain are both saved under the `data/calibration_data/from_calibration_chains` directory.
* `mcmc_acceptance_rate.csv` - the acceptance rate from the MCMC algorithm. Should be about 0.23 for the numbers of parameters (dimension) that we're dealing with here.

Note that calibrations of 20 million iterations will take multiple hours to complete.
* For BRICK on its own, this took about 8 hours on a standard desktop workstation (ca. 2020)
* For DOECLIM-BRICK and SNEASY-BRICK, this will take closer to 15 hours or so (using that same machine)

## Running the model hindcasts

This is done for the hindcast period 1850-2017 by using the `calibration/run_hindcast.jl` script, using `model_config=brick`, `doeclimbrick` and `sneasybrick` (three times). For the hindcast, no RCP scenario needs to be specified, because all of them follow historical radiative forcing/emissions trends up to 2005.

The standard set of parameters that are being used for the hindcast and projection simulations are the sub-sample of 10,000 from the MCMC calibration described above (`parameters_subsample.csv`). If you have a different parameter file that you want to run the hindcasts under, you will want to modify the section of `run_hindcast.jl` titled `Set paths for results files` (line 41).

This script will add to the date-sampled model configuration-specific directory that was constructed above (or came with the model codes). It will create a sub-directory called `hindcast_csv` which will be populated with CSV files that include the simulated hindcasts of the model output fields. Each of these names is appended with `model_config` (`brick`, `doeclimbrick`, or `sneasybrick`) and contains one hindcast simulation for each of the sets of parameters in the sub-sample for analysis. Rows correspond to different years (1850-2017 be default) and columns each correspond to different ensemble members.
* `hindcast_antarctic.csv` - contribution to global mean sea-level change from the Antarctic ice sheet (meters)
* `hindcast_greenland.csv` - contribution to global mean sea-level change from the Greenland ice sheet (meters)
* `hindcast_glaciers.csv` - contribution to global mean sea-level change from glaciers and ice caps (meters)
* `hindcast_landwater_storage_sl.csv` - contribution to sea-level change from land water storage (meters)
* `hindcast_gmsl.csv` - total global mean sea-level change (meters)
* `hindcast_ocean_heat.csv` - (DOECLIM- or SNEASY-BRICK only)
* `hindcast_temperature.csv` - (DOECLIM- or SNEASY-BRICK only)
* `hindcast_oceanco2.csv` - (SNEASY-BRICK only)
* `hindcast_co2.csv` - (SNEASY-BRICK only)
* `hindcast_MAP.csv` - all of the hindcast time series for the maximum _a posteriori_ set of parameters

## Running the model projections under different RCP scenarios

This is done for the period 1850-2300 (but can be modified to any period between 1765 and 2300) by using the `calibration/run_projections.jl` script, using `model_config=brick`, `doeclimbrick` or `sneasybrick` and `rcp_scenario="RCP26"`, `"RCP45"`, `"RCP60"`, or `"RCP85"`. Note that the RCP scenario forcing files are all the same until 2005.

This script will add to the date-sampled model configuration-specific directory that was constructed above (or came with the model codes). It will create a sub-directory called `projections_csv`, and a sub-directory within there that is specific to each RCP scenario used will be created. The projections files are analogous to the hindcast files that are generated, and will populate the `projections_csv/[RCP scenario]` directory.

## Creating forcing files for stand-alone BRICK

The forcing files for DOECLIM-BRICK (radiative forcing) and SNEASY-BRICK (emissions) are taken from the RCP database here (https://tntcat.iiasa.ac.at/RcpDb/dsd?Action=htmlpage&page=download) and the data repository of Malte Meinshausen here (http://www.pik-potsdam.de/~mmalte/rcps/data/).

For stand-alone BRICK, which requires temperature and ocean heat time series as forcing data, we use the time series for temperature and ocean heat uptake from the maximum _a posteriori_ simulations from the SNEASY-BRICK simulation ensembles described above. This is done in the `calibration/sneasy_make_brick_forcing.jl` script. This script creates the following files, where `xx` denotes the RCP scenario (`26`, `45`, `60`, or `85`), `yyyy` denotes the starting year of the forcing, and `YYYY` denotes the ending year of the forcing.
* `data/model_data/sneasy_oceanheat_RCPxx_yyyy_YYYY.csv`
* `data/model_data/sneasy_temperature_RCPxx_yyyy_YYYY.csv`

## Generating projections of local mean sea-level change

TODO

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
