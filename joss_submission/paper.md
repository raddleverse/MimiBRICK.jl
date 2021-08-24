---
title: 'MimiBRICK.jl: A Julia package for running and coupling the BRICK model for sea-level change in the Mimi integrated modeling framework'
tags:
  - Julia
  - sea level
  - climate
  - coastal
  - coupling
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
date: 24 August 2021
bibliography: paper.bib

---

# Summary

Assessment of policies to manage climate risks requires projections of future climate. 
For coastal risks, this includes projections of future global and local sea levels.
Major contributors to global sea-level change include glaciers and ice caps, thermal expansion, land water storage, and the Greenaland and Antarctic ice sheets.
**(TODO - finish general introduction to GMSL/LSL modeling, coastal hazards, and integrated modeling/Mimi idea)**

# Statement of need

The Building Blocks for Relevant Ice and Climate Knowledge (BRICK) model for sea-level rise [@Wong:2017] is a model for global and local mean sea-level change. 
Used in a number of recent assessments, including comparisons in the IPCC AR6 (TODO - citation). Relatedly, there is a need 

`Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

The original BRICK model was designed to be flexible and efficient **(TODO - continue)**... 
`MimiBRICK.jl` stays true to these design principles, and enhances the usability of the code by **(TODO - continue)**...

**(TODO - continue, talk about Mimi)**

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We gratefully acknowledge Corinne Hartin, Ryan Sriver, and Nathan Urban for their contributions.
**(TODO - add anyone else?)**

# References

