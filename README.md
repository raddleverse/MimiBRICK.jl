# MimiBRICK.jl

This repository will have the latest "standard" BRICK version and codes demonstrating its calibration in three configurations: (i) alone (forced by temperature), (ii) coupled to DOECLIM, and (iii) coupled to SNEASY. Standard calibration output will be provided, along with wrappers for Python, R, and Julia so that the model can be run in the same fashion from any of those three languages.

## Recommended Use

* If you would like to use previously published BRICK results as part of a new analysis, but do not necessarily want to re-run the model, then we recommend to go to the **Introduction-and-Library** repository. There, you will find links for each published BRICK study, including links for the calibrated model parameter data sets and the calibrated model projections for sea level, temperature, and any other relevant outputs.
* If you would like to run the model yourself, then you are in the right place!
  * Either fork the master branch from this repository or download the zipped file of codes.
  * Do your analysis
  * Please let us know of any new results that should be incorporated into the **Introduction-and-Library** repository! See the README.md file in that repository for more information.
  * Also please let us know of any model modifications and/or bug fixes that might usefully be incorporated into the main BRICK codes! Creating an "Issue" here is a great way to do that (top horizontal menu bar in the GitHub browser).

### License

GPL info to add

## How To Install Required Packages

This code was created using [Julia v1.5](https://julialang.org/downloads/) and requires several Julia packages.

(1) To install these packages, first enter the package manager by hitting the `]` key in the Julia console. Once in the package manager, run the following code:

```julia
add CSVFiles  
add DataFrames  
add Distributions
add KernelDensity
add Mimi  
add NetCDF
add RobustAdaptiveMetropolisSampler
```
(2) While still in the package manager, run the following line to install the Mimi implementation of SNEASY:

```julia
add https://github.com/anthofflab/MimiSNEASY.jl.git
```

(3) Run the following line to install the Mimi implementation of BRICK:

```julia
add https://github.com/FrankErrickson/MimiBRICK.jl.git
```

(4) To exit back to Julia, hit the `backspace` key.

## Run the Baseline SNEASY+BRICK Calibration

**(here is where we may deviate from what Frank had previously)**

(1) First, [Clone or download](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository) the `brick_scc_paper` Git repository. Once this is on your computer, set this folder as your working directory in the Julia console.

(2) In lines 24-41 of the `scr/calibrate_sneasybrick_for_tony.jl` file, set the folder name to save results and the number of samples to take during the MCMC calibration. Save this file.

(3) Run the following line of code to carry out the SNEASY+BRICK calibration and automatically save results.

```julia
include("scr/calibrate_sneasybrick_for_tony.jl")
```

## Description of SNEASY+BRICK Baseline Calibration Files

(1) `calibration/calibration_helper_functions.jl`: Contains various functions that are useful for the model calibration.

(2) `calibration/run_historic_models/run_sneasy_brick_historic_climate.jl`: Creates an instance of SNEASY+BRICK that will automatically update model projections over the hindcast period when passing in a new set of parameter values (mostly there to make calibration code run faster).

(3) `calibration/create_log_posterior_sneasy_brick.jl`: Creates functions to calculate the prior, likelihood, and posterior for SNEASY+BRICK that can then be passed into the MCMC calibration.

(4) `src/calibrate_sneasybrick_for_tony.jl`: Contains a few model settings at the top (length of MCMC chain, name of folder to save results, etc.), and then otherwise will load all of the necessary files and carries out the SNEASY+BRICK calibration.
