![](https://github.com/raddleverse/MimiBRICK.jl/workflows/Run%20CI%20on%20master/badge.svg)
[![codecov](https://codecov.io/gh/raddleverse/MimiBRICK.jl/branch/master/graph/badge.svg?token=H7SJB47W5V)](https://codecov.io/gh/raddleverse/MimiBRICK.jl)

# MimiBRICK.jl

This repository will have the latest "standard" BRICK version and codes demonstrating its calibration in three configurations: (i) alone (forced by temperature and ocean heat), (ii) coupled to DOECLIM, and (iii) coupled to SNEASY. Standard calibration output will be provided as well as examples and scripts for downscaling the projections to a 1-degree grid for local sea-level rise.

<br>

## Statement of Need

MimiBRICK.jl is an implementation of the [Building Blocks for Relevant Ice and Climate Knowledge (BRICK) semi-empirical model for sea-level change](https://doi.org/10.5194/gmd-10-2741-2017) in the Mimi integrated modeling framework (https://www.mimiframework.org/). The Mimi modeling framework is a coding platform that facilitates coupling models and running coupled modeling experiments. MimiBRICK.jl is flexible, efficient, and modular, to facilitate incorporating BRICK into coupled models and integrated assessments of climate impacts in a modular fashion to provide global average as well as local sea-level projections. This focus on tight model coupling and integrated modeling is a key feature of MimiBRICK.jl and broader Mimi modeling framework.

This implementation includes examples for using observational data to calibrate the model, as well as various configurations in which MimiBRICK.jl is coupled to other climate model components. For users who do not wish to re-run computationally intensive model calibration algorithms, this implementation also includes scripts for using existing calibration output for standard future climate change scenarios, and examples downscaling these global projections for assessments of local impacts. Pre-run model calibration and simulation output can be found in the [accompanying Zenodo repository](https://zenodo.org/record/6626335).

<br>

## Recommended Use

* If you would like to use previously published BRICK results as part of a new analysis, but do not necessarily want to re-run the model, then we recommend to go to the [**Introduction-and-Library** repository](https://github.com/raddleverse/Introduction-and-Library). There, you will find links for each published BRICK study, including links for the calibrated model parameter data sets and the calibrated model projections for sea level, temperature, and any other relevant outputs. Accompanying this package, in the `examples` directory, you will find a Zenodo repository linked. This repository has a set of calibration output and model ensembles that can be used without re-calibrating the model.
* If you would like to run the model yourself, then you are in the right place!
  * It is assumed that users will not clone/download this Git repository. Instead, you can add and use the package as described below. Those commands and the commands to set up and run the model can be executed from the directories that contain the rest of your project codes.
  * Load the MimiBRICK package (or if you wish to edit the package either fork the master branch from this repository or download the zipped file of codes)
  * Do your analysis
  * Please let us know of any new results that should be incorporated into the [Introduction-and-Library repository](https://github.com/raddleverse/Introduction-and-Library)! See the README.md file in that repository for more information.
  * Also please let us know of any model modifications and/or bug fixes that might usefully be incorporated into the main BRICK codes. Creating an "Issue" here is a great way to do that (top horizontal menu bar in the GitHub browser).

<br>

## Installation and Examples

Installation instructions and example use cases can be found in the GitHub Pages site: https://raddleverse.github.io/MimiBRICK.jl/

<br>

## Author Contributions

* TW: initial model development, software development, model calibration and validation, conceptualization, projection direction and overall management
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
