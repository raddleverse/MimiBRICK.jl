---
title: 'MimiBRICK.jl: A Julia package for the BRICK model for sea-level change in the Mimi integrated modeling framework'
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
 - name: Thayer School of Engineering, Dartmouth College, USA
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
`MimiBRICK.jl` is an implementation of the Building Blocks for Relevant Ice and Climate Knowledge (BRICK) semi-empirical model for sea-level change [@Wong2017_brick] in the Mimi framework.
`MimiBRICK.jl` is flexible, efficient, and modular, to facilitate incorporating BRICK into coupled models and integrated assessments of climate impacts in a modular fashion [@national_academies_of_sciences_valuing_2017].
This implementation includes examples for using observational data to calibrate the model, as well as various configurations in which `MimiBRICK` is coupled to other climate model components.
For users who do not wish to re-run computationally intensive model calibration algorithms, this implementation also includes code for using existing calibration output for standard future climate change scenarios, and examples downscaling these global projections for assessments of local impacts.

# Summary

The BRICK semi-empirical model for sea-level rise [@Wong2017_brick] is a model for global and local mean sea-level change.
The core model includes component sub-models for the major contributors to global mean sea-level change (glaciers and ice caps, thermal expansion, land water storage, and the Greenland and Antarctic ice sheets).
The resulting global mean sea levels can be downscaled via a data set that represents the "fingerprint" of each sea-level component on local mean sea level [@slangen2014].
In this way, BRICK provides useful information about local sea-level changes, including characterizing uncertainties and being flexible and efficient enough to resolve high-risk upper tails of probability distributions.
BRICK has been used in a number of recent assessments, including for examining the impacts of sea-level rise as a constraint on estimates of climate sensitivity [@Vega-Westhoff2018], estimates of deep uncertainty in coastal flood risk [@Ruckert2019], and most recently was included in comparisons of sea-level projections in the Sixth Assessment Report of the Intergovernmental Panel on Climate Change [@ipccar6ch9].

By working with annual global mean temperatures and sea levels, BRICK is suitable for embedding within and coupling to other models for climate change and its impacts. `MimiBRICK.jl`, written in compliance with the Mimi integrated modeling framework, facilitates incorporating BRICK into larger-scale coupled modeling efforts. The `MimiBRICK.jl` repository includes three such examples: (i) standalone BRICK, which takes as input temperature and ocean heat uptake; (ii) BRICK coupled to a simple one-dimensional Diffusion-Ocean-Energy Climate model (DOECLIM) coupled to BRICK, which takes as input radiative forcing scenarios such as the standard Representative Concentration Pathway (RCP) scenarios; and (iii) BRICK coupled to a Simple Nonlinear Earth System model (SNEASY) coupled to BRICK, which takes as input radiative forcing and greenhouse gas emissions and concentration scenarios (such as the RCP scenarios).

The standalone BRICK model requires as input annual mean time series for global mean surface temperature and ocean heat uptake.  In the DOECLIM-BRICK and SNEASY-BRICK configurations, those temperatures and ocean heat uptake inputs are provided to BRICK through output from the DOECLIM and SNEASY models. Model output and parameter values can be explored using the Mimi `explore()` function (Figure 1). The `explore()` function allows users to easily view and zoom in on different features in the model simulation set up or results. Being coded in the Mimi style enables the user to more easily couple BRICK to the suite of other models already implemented in Mimi, and builds on the extensive documentation and community support online (https://www.mimiframework.org/Mimi.jl/stable/).

# Acknowledgements

We gratefully acknowledge Corinne Hartin, Ryan Sriver, and Nathan Urban for their contributions.
**(TODO - add anyone else?)**

# References
