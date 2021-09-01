# MimiBRICK.jl

This repository will have the latest "standard" BRICK version and codes demonstrating its calibration in three configurations: (i) alone (forced by temperature), (ii) coupled to DOECLIM, and (iii) coupled to SNEASY. Standard calibration output will be provided as well as examples and scripts for downscaling the projections to a 1-degree grid for local sea-level rise.

## Recommended Use

* If you would like to use previously published BRICK results as part of a new analysis, but do not necessarily want to re-run the model, then we recommend to go to the **Introduction-and-Library** repository. There, you will find links for each published BRICK study, including links for the calibrated model parameter data sets and the calibrated model projections for sea level, temperature, and any other relevant outputs.
* If you would like to run the model yourself, then you are in the right place!
  * Either fork the master branch from this repository or download the zipped file of codes.
  * Do your analysis
  * Please let us know of any new results that should be incorporated into the **Introduction-and-Library** repository! See the README.md file in that repository for more information.
  * Also please let us know of any model modifications and/or bug fixes that might usefully be incorporated into the main BRICK codes. Creating an "Issue" here is a great way to do that (top horizontal menu bar in the GitHub browser).

### License

Copyright 2021 Tony Wong, (todo)

This file is part of MimiBRICK.jl (Building blocks for Relevant Ice and Climate Knowledge). MimiBRICK.jl is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

MimiBRICK.jl is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with MimiBRICK.jl (`LICENSE.md)`). If not, see http://www.gnu.org/licenses/.


## How To Install Required Packages

This code was created using [Julia v1.6](https://julialang.org/downloads/) and requires several Julia packages.

(1) To install these packages, first enter the package manager by hitting the `]` key in the Julia console. Once in the package manager, run the following code:

```julia
add CSVFiles
add DataFrames
add Distributions
add KernelDensity
add Mimi
add MimiSNEASY
add NetCDF
add RobustAdaptiveMetropolisSampler
```
(2) While still in the package manager, run the following line to install the Mimi implementation of SNEASY:

```julia
add https://github.com/anthofflab/MimiSNEASY.jl.git
```

(3) Run the following line to install the Mimi implementation of BRICK:

```julia
add https://github.com/BRICK-SLR/MimiBRICK.jl.git
```

(4) To exit back to Julia, hit the `backspace` key.

## Running baseline cases with default parameters

As a preliminary step, you should activate the Julia environment by doing the following in a Julia console:
1. Navigate to the

### BRICK standalone (with temperature and ocean heat uptake exogenous forcing)

This is the first test that is done in `test/runtests.jl`. Since it does not require DOECLIM or SNEASY, you can run BRICK using temperature and ocean heat uptake forcing data by running in the Julia console:
```
using MimiBRICK
m = MimiBRICK.get_model()
run(m)
```

You can plot the output fields in the model object `m` using (for example) the `Plots` Julia plotting package. First, let's grab the years over which the model was run. This is a dimension in the model. We can retrieve it by using the `dim_keys()` function, from the `Mimi` package.
```
using Mimi
years = dim_keys(m, :time)
```

Then we can load the `Plots` package and make a figure of the global mean sea-level change. Note that the first argument into the `m` object specifies the component of our model, and the second argument specifies the field. Here, we are grabbing the `sea_level_rise` field from the `global_sea_level` component.
```
using Plots
plot(years, m[:global_sea_level, :sea_level_rise])
```

Mimi also offers an explorer window to check these model output fields out. To use this, we need to load the `Mimi` package (if you haven't already).
```
using Mimi
```

Then, we can open the explorer.
```
explore(m)
```

This should open a window labeled "Mimi Explorer Window". On the left, there should be two vertically-stacked boxes. The top box is labeled "Components" and the bottom box is labeled "Data". To view some of the model output from our out-of-box BRICK simulation, you must first pick one of the Components from the top box, and then a Variable field out of the bottom box. For example, to view the global mean sea level model projections, select `global_sea_level` from the Components box, and `global_sea_level : sea_level_rise` from the Variables box. A plot of this variable should appear in the right column of plot boxes. The top box is static; the bottom box you can interact with to zoom in on different portions of the figure.

More information about exploring Mimi model results can be found in the [Mimi Framework How-To guides online](https://www.mimiframework.org/Mimi.jl/stable/howto/howto_2/).

### BRICK+DOECLIM (with radiative forcing)

TODO

### BRICK+SNEASY (with emissions forcing)

TODO

## Running baseline cases with parameters provided via CSV file

TODO

## Running simulations with parameters from a previous ensemble

TODO

## Run the Baseline SNEASY+BRICK Calibration

**(TODO - modify descriptions, what Frank had previously)**

(1) First, [Clone or download](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository) the `brick_scc_paper` Git repository. Once this is on your computer, set this folder as your working directory in the Julia console.

(2) In lines 24-41 of the `scr/calibrate_sneasybrick_for_tony.jl` file, set the folder name to save results and the number of samples to take during the MCMC calibration. Save this file.

(3) Run the following line of code to carry out the SNEASY+BRICK calibration and automatically save results.

```julia
include("scr/calibrate_sneasybrick_for_tony.jl")
```

## Description of SNEASY+BRICK Baseline Calibration Files

**(TODO - modify descriptions, what Frank had previously)**

(1) `calibration/calibration_helper_functions.jl`: Contains various functions that are useful for the model calibration.

(2) `calibration/run_historic_models/run_sneasy_brick_historic_climate.jl`: Creates an instance of SNEASY+BRICK that will automatically update model projections over the hindcast period when passing in a new set of parameter values (mostly there to make calibration code run faster).

(3) `calibration/create_log_posterior_sneasy_brick.jl`: Creates functions to calculate the prior, likelihood, and posterior for SNEASY+BRICK that can then be passed into the MCMC calibration.

(4) `src/calibrate_sneasybrick_for_tony.jl`: Contains a few model settings at the top (length of MCMC chain, name of folder to save results, etc.), and then otherwise will load all of the necessary files and carries out the SNEASY+BRICK calibration.

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
