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
    affiliation: "5,7"
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
 - name: Department of Hydraulic Engineering, Faculty of Civil Engineering and Geosciences, Delft University of Technology, The Netherlands
   index: 7
date: 17 June 2022
bibliography: paper.bib

---

# Statement of need

Assessing strategies to manage climate risks requires sound projections of future climate. For coastal risks, this includes projections of future sea levels. Major contributors to sea-level change include glaciers and ice caps, thermal expansion, land water storage, and the Greenland and Antarctic ice sheets. Characterizing coastal hazards and managing the associated risks requires resolving the tails of distributions. Semi-empirical models for sea-level rise offer a computationally efficient method for characterizing uncertainties in future coastal hazards and fusing observational data with models [@Kopp2017; @Mengel2016; @Nauels2017; @Wong2017_brick]. These models can also be flexible and modular, which enables their use in integrated frameworks for assessing climate damages and examining the efficacy of climate risk management policies. The Mimi integrated modeling framework (https://www.mimiframework.org/) is a coding platform that facilitates coupling models and running coupled modeling experiments.`MimiBRICK.jl` is an implementation of the Building Blocks for Relevant Ice and Climate Knowledge (BRICK) semi-empirical model for sea-level change [@Wong2017_brick] in the Mimi framework. `MimiBRICK.jl` is flexible, efficient, and modular, to facilitate incorporating BRICK into coupled models and integrated assessments of climate impacts in a modular fashion to provide global average as well as local sea-level projections [@national_academies_of_sciences_valuing_2017].

This implementation includes examples for using observational data to calibrate the model, as well as various configurations in which `MimiBRICK` is coupled to other climate model components. For users who do not wish to re-run computationally intensive model calibration algorithms, this implementation also includes scripts for using existing calibration output for standard future climate change scenarios, and examples downscaling these global projections for assessments of local impacts.

# Summary

BRICK is a semi-empirical model for global and local mean sea-level change [@Wong2017_brick]. The core model includes component sub-models for the major contributors to global mean sea-level change - glaciers and ice caps, thermal expansion, land water storage, and the Greenland and Antarctic ice sheets. The resulting global mean sea levels can be downscaled via a data set that represents the "fingerprint" of each sea-level component on local mean sea level [@slangen2014]. In this way, BRICK provides information about local sea-level changes, including characterizations of key uncertainties. BRICK is flexible and efficient enough to resolve high-risk upper tails of probability distributions. BRICK has been used in a number of recent assessments, including for examining the impacts of sea-level rise as a constraint on estimates of climate sensitivity [@Vega-Westhoff2018], estimates of deep uncertainty in coastal flood risk [@Ruckert2019], and most recently was included in comparisons of sea-level projections in the Sixth Assessment Report of the Intergovernmental Panel on Climate Change [@ipccar6ch9].

By working with annual global mean temperatures and sea levels, BRICK is suitable for embedding within, and coupling to, other models for climate change and its impacts. `MimiBRICK.jl` is written in compliance with the Mimi integrated modeling framework to facilitate incorporating BRICK into larger-scale coupled modeling efforts. The `MimiBRICK.jl` repository includes three such examples: (i) standalone BRICK, which takes as input temperature and ocean heat uptake; (ii) BRICK coupled to a simple one-dimensional Diffusion-Ocean-Energy Climate model (DOECLIM), which takes as input radiative forcing scenarios such as the standard Representative Concentration Pathway (RCP) scenarios; and (iii) BRICK coupled to a Simple Nonlinear Earth System model (SNEASY), which takes as input radiative forcing and greenhouse gas emissions and concentration scenarios (such as the RCP scenarios).

The standalone BRICK model requires as input annual mean time series for global mean surface temperature and oceanic heat uptake. In the DOECLIM-BRICK and SNEASY-BRICK configurations, those temperatures and ocean heat uptake inputs are provided to BRICK through output from the DOECLIM and SNEASY models. In the examples provided in the `MimiBRICK.jl` repository, the temperature time series is from the maximum _a posteriori_ simulation from the ensemble run using the SNEASY-BRICK configuration. In all Mimi modeling experiments, model output and parameter values can be explored using the Mimi `explore()` function (\autoref{fig:explorer}). The `explore()` function allows users to easily view and zoom in on different features in the model simulation set up or results. Being coded in the Mimi style enables the user to more easily couple BRICK to the suite of other models already implemented in Mimi, and builds on the extensive documentation and community support online (https://www.mimiframework.org/Mimi.jl/stable/).

![The Mimi `explore()` function allows users to interactively explore the coupled model output variables and parameters associated with each model component.\label{fig:explorer}](mimi_explorer.png)

# Acknowledgements

We gratefully acknowledge Corinne Hartin, Matthew Hoffman, Radley Powers, Ryan Sriver, Nathan Urban, and Benjamin Vega-Westhoff for valuable contributions. This work was co-supported in part by the Penn State Center for Climate Risk Management and the Thayer School of Engineering. All errors and opinions are those of the authors and not of the supporting entities.

Software License: The MimiBRICK.jl code is distributed under GNU general public license. The authors do not assume responsibility for any (mis)use of the provided code.

# References
