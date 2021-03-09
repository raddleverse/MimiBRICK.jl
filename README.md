# MimiBRICK.jl

### How To Install Required Packages

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

### Run the Baseline SNEASY+BRICK Calibration

**(here is where we may deviate from what Frank had previously)**

(1) First, [Clone or download](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository) the `brick_scc_paper` Git repository. Once this is on your computer, set this folder as your working directory in the Julia console.

(2) In lines 24-41 of the `scr/calibrate_sneasybrick_for_tony.jl` file, set the folder name to save results and the number of samples to take during the MCMC calibration. Save this file.

(3) Run the following line of code to carry out the SNEASY+BRICK calibration and automatically save results.

```julia
include("scr/calibrate_sneasybrick_for_tony.jl")
```

### Description of SNEASY+BRICK Baseline Calibration Files

(1) `calibration/calibration_helper_functions.jl`: Contains various functions that are useful for the model calibration. 

(2) `calibration/run_historic_models/run_sneasy_brick_historic_climate.jl`: Creates an instance of SNEASY+BRICK that will automatically update model projections over the hindcast period when passing in a new set of parameter values (mostly there to make calibration code run faster).

(3) `calibration/create_log_posterior_sneasy_brick.jl`: Creates functions to calculate the prior, likelihood, and posterior for SNEASY+BRICK that can then be passed into the MCMC calibration.

(4) `src/calibrate_sneasybrick_for_tony.jl`: Contains a few model settings at the top (length of MCMC chain, name of folder to save results, etc.), and then otherwise will load all of the necessary files and carries out the SNEASY+BRICK calibration.
