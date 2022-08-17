![](https://github.com/raddleverse/MimiBRICK.jl/workflows/Run%20CI%20on%20master/badge.svg)
[![codecov](https://codecov.io/gh/raddleverse/MimiBRICK.jl/branch/master/graph/badge.svg?token=H7SJB47W5V)](https://codecov.io/gh/raddleverse/MimiBRICK.jl)

# MimiBRICK.jl

This repository will have the latest "standard" BRICK version and codes demonstrating its calibration in three configurations: (i) alone (forced by temperature and ocean heat), (ii) coupled to DOECLIM, and (iii) coupled to SNEASY. Standard calibration output will be provided as well as examples and scripts for downscaling the projections to a 1-degree grid for local sea-level rise.

<br>

## Statement of Need

MimiBRICK.jl is an implementation of the [Building Blocks for Relevant Ice and Climate Knowledge (BRICK) semi-empirical model for sea-level change](https://doi.org/10.5194/gmd-10-2741-2017) in the Mimi integrated modeling framework (https://www.mimiframework.org/). The Mimi modeling framework is a coding platform that facilitates coupling models and running coupled modeling experiments. MimiBRICK.jl is flexible, efficient, and modular, to facilitate incorporating BRICK into coupled models and integrated assessments of climate impacts in a modular fashion to provide global average as well as local sea-level projections. This focus on tight model coupling and integrated modeling is a key feature of MimiBRICK.jl and broader Mimi modeling framework.

This implementation includes examples for using observational data to calibrate the model, as well as various configurations in which MimiBRICK.jl is coupled to other climate model components. For users who do not wish to re-run computationally intensive model calibration algorithms, this implementation also includes scripts for using existing calibration output for standard future climate change scenarios, and examples downscaling these global projections for assessments of local impacts. Pre-run model calibration and simulation output can be found in the [accompanying Zenodo repository](https://zenodo.org/record/6626335).

<br>

## Structure

* `/calibration` - functions for setting up model forcing for calibration
* `/data` - calibration and forcing scenario data
* `/docs` - files for GitHub Pages site documentation
* `/examples` - Jupyter notebooks demonstrating the workflows for model calibration, downscaling global sea-level change to local, and making model hindcasts and projections; we anticipate that it will be useful to take these notebooks and modifying them to fit users' own use cases
* `/joss_submission` - files associated with the MimiBRICK.jl Journal of Open Source Software submission
* `/src` - functions for the actual component submodels of BRICK, for configuring these models as a combined coupled BRICK model, for performing the downscaling to local sea-level change, and for performing the model calibration; includes the likelihood function configuration
* `/test` - contains tests for out-of-box model configurations, testing a small model calibration, and downscaling to local sea-level change; used for continuous integration testing

### Substituting/Modifying Component Sub-models

Below is a summary of the changes that would be needed to modify a component sub-model of MimiBRICK. A much more comprehensive discussion of running, modifying old, and creating new models within the Mimi coupled modeling framework can be found at the [Mimi Documentation](https://www.mimiframework.org/Mimi.jl/stable/).

* Create/modify the source code for the new/modified component in `/src/components`. See the Design Pattern section below for a stencil. Pattern-matching with the existing component sub-models is encouraged!
* In `/src/MimiBRICK.jl`…
    * Add an `include` statement if you have created a new file for your model component
    * Add/modify `update_param` calls, component name, and parameter names/default values
    * Connect parameters across components using `connect_param` and `add_shared_param` near the bottom of the `brick` module, as needed
* In `/src/calibration/run_historic_models/run_brick_historic_climate.jl`…
    * Add/modify the parameter names and `update_param` calls to match what is in the source code for your modified component and in `/src/MimiBRICK.jl`
* In `/src/calibration/run_hindcast.jl` and `/src/calibration/run_projections.jl`…
    * Modify the `ar1_noise_xxx` and `obs_error_xxx` variables, as needed, depending on whether the modified component uses a similar or new residual model. Might also require modifying at “Statistical noise models” if a new residual model is used
    * Modify `update_param` calls associated with the modified component (similarly to in `/src/MimiBRICK.jl`)
    * If your modification create new output that you would like to write to files, you may want to modify at “Save output” as well to include any new fields that you want to write
* If your altered the model components and the naming conventions for the model output CSV files, then you may also need to modify the `/src/downscale.jl` file to generate estimates of local sea-level change.
* If you would like to modify the calibration data used, then you should modify…
    * the function `load_calibration_data` within `/src/calibration/calibration_helper_functions.jl`
    * the appropriate `create_log_posterior_xxx.jl` script within `/src/calibration/create_log_posteriors`

### Design Pattern

Typical Mimi component models use the pattern depicted below. Most notably, models are run timestep-by-timestep using a `run_timestep` function that is defined for each model component. `run_timestep` has four input arguments:
* `p` - model parameters (see [Mimi Documentation on parameters and variables](https://www.mimiframework.org/Mimi.jl/stable/howto/howto_5/) for more information)
* `v` - model variables
* `d` - model dimensions (e.g., `d.regions` or `d.time`)
* `t` - timestep (see [Mimi Documentation on the timestep types](https://www.mimiframework.org/Mimi.jl/stable/howto/howto_4/) for more information)

`run_timestep` typically does not have a return value. Rather, this function modifies the shared model variables in `v`. For example, within `MimiBRICK.jl`, sea levels for the glaciers and ice caps component are computed and stored in `v.gsic_sea_level` within the glaciers and ice caps `run_timestep` function.

More information about defining new model components in Mimi can be found in the [Mimi Documentation](https://www.mimiframework.org/Mimi.jl/stable/tutorials/tutorial_4/).

```
using Mimi

@defcomp component_name begin

    # --------------------
    # Model Parameters
    # --------------------

    compshortname_scalarparam = Parameter()                   # description (units)
    compshortname_timeseriesparam = Parameter(index=[time])   # description (units)

    # --------------------
    # Model Variables
    # --------------------

    compshortname_timeseriesvariable = Variable(index=[time]) # description (units)

    # --------------------
    # Model Equations
    # --------------------
    
    function run_timestep(p, v, d, t)
    
        # what does a single time step within the sub-model component do?
      
    end
    
end
```

### Recommended Use

* If you would like to use previously published BRICK results as part of a new analysis, but do not necessarily want to re-run the model, then we recommend to go to the [accompanying Zenodo repository](https://zenodo.org/record/6626335). In the near future, we will compile a library of other previously published studies using BRICK, including any other relevant outputs. 
* If you would like to run the model yourself, then you are in the right place!
  * It is assumed that users will not clone/download this Git repository. Instead, you can add and use the package as described below. Those commands and the commands to set up and run the model can be executed from the directories that contain the rest of your project codes.
  * Load the MimiBRICK package (or if you wish to edit the package either fork the master branch from this repository or download the zipped file of codes)
  * Do your analysis
  * Also please let us know of any model modifications and/or bug fixes that might usefully be incorporated into the main BRICK codes. Creating an "Issue" here is a great way to do that (top horizontal menu bar in the GitHub browser).

### Installation and Examples

This code was created using Julia v1.6 and requires several Julia packages. It is recommended that you use Julia v1.6 (or later). Julia may be downloaded from http://julialang.org/downloads/.

(1) Run the following line to install the Mimi implementation of BRICK:

```
]
add MimiBRICK
```

(2) To exit back to Julia, hit the backspace key.

Further instructions and example use cases can be found in the `examples` directory and on the GitHub Pages site https://raddleverse.github.io/MimiBRICK.jl/.

<br>

## Contributions

Contributions to Mimi codes are most welcome! More information can be found at in the Mimi Framework site at https://www.mimiframework.org/. Users are encouraged to engage with the development team on GitHub and in the [Mimi Framework Forum](https://forum.mimiframework.org/).

### Contributions to the initial v1.0 of MimiBRICK.jl

* TW: initial model development, software development, model calibration and validation, conceptualization, projection direction, and overall management
* LR: software development, package maintenance, conceptualization
* FE: software development, model calibration and validation, conceptualization
* VS: software testing, model calibration and validation, conceptualization
* AB: initial model development, software testing, conceptualization
* KK: software testing, conceptualization
* DA: software development, package maintenance, conceptualization, project direction

<br>

## License

Copyright 2022 Tony Wong, Lisa Rennels, Frank Errickson, Vivek Srikrishnan, Alexander Bakker, Klaus Keller, and David Anthoff

This file and codes in this repository are part of MimiBRICK.jl (Building blocks for Relevant Ice and Climate Knowledge). MimiBRICK.jl is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

MimiBRICK.jl is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with MimiBRICK.jl (`LICENSE.md)`). If not, see http://www.gnu.org/licenses/.
