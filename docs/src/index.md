# MimiBRICK.jl

MimiBRICK.jl is an implementation of the [Building Blocks for Relevant Ice and Climate Knowledge (BRICK) semi-empirical model for sea-level change](https://doi.org/10.5194/gmd-10-2741-2017) in the Mimi integrated modeling framework (https://www.mimiframework.org/). The Mimi modeling framework is a coding platform that facilitates coupling models and running coupled modeling experiments. MimiBRICK.jl is flexible, efficient, and modular, to facilitate incorporating BRICK into coupled models and integrated assessments of climate impacts in a modular fashion to provide global average as well as local sea-level projections. This focus on tight model coupling and integrated modeling is a key feature of MimiBRICK.jl and broader Mimi modeling framework.

This implementation includes examples for using observational data to calibrate the model, as well as various configurations in which MimiBRICK.jl is coupled to other climate model components. For users who do not wish to re-run computationally intensive model calibration algorithms, this implementation also includes scripts for using existing calibration output for standard future climate change scenarios, and examples downscaling these global projections for assessments of local impacts. Pre-run model calibration and simulation output can be found in the [accompanying Zenodo repository](https://zenodo.org/record/6626335).

```@autodocs
Modules = [MimiBRICK]
```



## Author Contributions

* TW: initial model development, software development, model calibration and validation, conceptualization, projection direction and overall management
* LR: software development, package maintenance, conceptualization
* FE: software development, model calibration and validation, conceptualization
* VS: software testing, model calibration and validation, conceptualization
* AB: initial model development, software testing, conceptualization
* KK: software testing, conceptualization
* DA: software development, package maintenance, conceptualization, project direction



## License

Copyright 2022 Tony Wong, Lisa Rennels, Frank Errickson, Vivek Srikrishnan, Alexander Bakker, Klaus Keller, and David Anthoff

This file and codes in this repository are part of MimiBRICK.jl (Building blocks for Relevant Ice and Climate Knowledge). MimiBRICK.jl is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

MimiBRICK.jl is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with MimiBRICK.jl (`LICENSE.md)`). If not, see http://www.gnu.org/licenses/.
