---
title: 'MimiBRICK.jl: A Julia package for running and coupling the BRICK model for sea-level change in the Mimi integrated modeling framework'
tags:
  - sea level
  - climate
  - coastal
  - model coupling
  - Julia
authors:
  - name: Tony E. Wong
    orcid: 0000-0002-7304-3883
    affiliation: 1
  - name: Lisa Rennels
    orcid: 0000-0002-5307-9336
    affiliation: 2
  - name: Frank Errickson
    orcid: 0000-0003-2195-9424
    affiliation: 3
  - name: Vivek Srikrishnan
    orcid: 0000-0003-0049-3805
    affiliation: 4
  - name: Alexander Bakker
    orcid: 0000-0002-1017-7665
    affiliation: 5
  - name: Klaus Keller
    orcid: 0000-0002-5451-8687
    affiliation: 6
  - name: David Anthoff
    orcid: 0000-0001-9319-2109
    affiliation: 2
affiliations:
 - name: School of Mathematical Sciences, Rochester Institute of Technology, USA
   index: 1
 - name: Energy and Resources Group, University of California, Berkeley, USA
   index: 2
 - name: School of Public and International Affairs, Princeton University, USA
   index: 3
 - name: Department of Biological and Environmental Engineering, Cornell University, USA
   index: 4
 - name: Rijkswaterstaat, Ministry of Infrastructure and Water Management, The Netherlands
   index: 5
 - name: Department of Geosciences, Pennsylvania State University, USA
   index: 6
date: 27 August 2021
bibliography: paper.bib

---

# Statement of need

Assessment of policies to manage climate risks requires projections of future climate.
For coastal risks, this includes projections of future global and local sea levels.
Major contributors to global sea-level change include glaciers and ice caps, thermal expansion, land water storage, and the Greenland and Antarctic ice sheets.
Characterizing coastal hazards and managing the associated risks requires resolving the tails of distributions.
Semi-empirical models for sea-level rise offer a computationally efficient method for characterizing uncertainties in future coastal hazards and fusing observational data with models [@Kopp2017; @Mengel2016; @Nauels2017; @Wong2017_brick].
These models can also be flexible and modular, which enables their use in integrated frameworks for assessing climate damages and examining the efficacy of climate risk management policies.
The Mimi integrated modeling framework (https://www.mimiframework.org/) is a coding platform that facilitates coupling models and running coupled modeling experiments.
`MimiBRICK.jl` is an implementation of the Building Blocks for Relevant Ice and Climate Knowledge (BRICK) semi-empirical model for sea-level change [@Wong2017_brick] in the Mimi modeling framework.
`MimiBRICK.jl` is flexible and efficient, purposefully structured to be coupled into integrated assessments of climate impacts.
This implementation includes examples for using observational data to calibrate the model, as well as various configurations in which `MimiBRICK` is coupled to other climate model components.
For users who do not wish to re-run computationally intensive model calibration algorithms, this implementation also includes code for using existing calibration output for standard future climate change scenarios, and examples downscaling these global projections for assessments of local impacts.

# Summary

The BRICK semi-empirical model for sea-level rise [@Wong2017_brick] is a model for global and local mean sea-level change.
Used in a number of recent assessments, including comparisons in the IPCC AR6 (TODO - citation).

The original BRICK model was designed to be flexible and efficient **(TODO - continue)**...
`MimiBRICK.jl` stays true to these design principles, and enhances the usability of the code by **(TODO - continue)**...

`MimiBRICK.jl` is an implementation of the BRICK model in the Mimi integrated modeling framework

**(TODO - continue, talk about Mimi)**

# Acknowledgements

We gratefully acknowledge Corinne Hartin, Ryan Sriver, and Nathan Urban for their contributions.
**(TODO - add anyone else?)**

# References
