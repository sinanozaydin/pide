---
title: 'pide: Petrophysical Interpretation tools for geoDynamic Exploration.'
tags:
  - Python3
  - magnetotellurics
  - petrophysics
  - electrical conductivity
  - seismic velocity
  - experimental petrology
  - geodynamics
  - geophysics
  - synthetic modelling
authors:
  - name: Sinan Ozaydin
    orcid: 0000-0002-4532-9980
    equal-contrib: true
    corresponding: true 
    affiliation: 1 
  - name: Lu Li
    orcid: 0000-0001-6165-2828
    affiliation: 2
  - name: Utpal Singh
    orcid: 0000-0001-8304-5615
    affiliation: 1
  - name: Patrice F. Rey
    orcid: 0000-0002-1767-8593
    affiliation: 1
  - name: Maria Constanza Manassero
    orcid: 0000-0002-3178-6523
    affiliation: 3

affiliations:
 - name: School of Geosciences, University of Sydney, Sydney, Australia.
   index: 1
 - name: School of Earth Sciences, University of Western Australia, Perth, Australia.
   index: 2
 - name: School of Natural Sciences (Physics), University of Tasmania, Hobart, Australia.
   index: 3

date: 6 May 2024
bibliography: paper.bib
---

# Summary

`pide` is a Python3 library for calculating geophysical parameters (e.g., electrical conductivity, seismic velocity), employing the results from experimental petrology, mineral/rock physics, and thermomechanical modelling studies. `pide` can calculate the theoretical electrical conductivity of any earth material that exists in the literature. `pide` can also calculate seismic velocity utilising the external 'sister' library `SAnTex`. Using these theoretical calculations, users can utilise inversion modules to decode geophysical anomalies compositionally or convert thermomechanical models into geophysical observables. With a given spatial mapping of earth materials, which can preferentially be loaded from a thermomechanical model, `pide`  can be used to build synthetic electrical conductivity and seismic velocity models and generate gravity and magnetic anomalies. Moreover, it is built as a modular tool, so users can easily build their functions.

# Statement of need
Given the inherent heterogeneity and complexity of Earth systems, geophysical tomographies often yield complex 2D and 3D images that are challenging to interpret. To enhance their interpretations, researchers commonly turn to experimental petrology and mineral physics, covering various geophysical properties, including electrical conductivity, magnetic susceptibility, seismic velocity, and rheology. These properties are sensitive to phase transitions, partial melting, major and trace elements partitioning, mineral solubilities, and phase-mixing models. Numerous specialized tools have been designed to address specific properties, many of which feature graphical interfaces [@Ozaydin2020; @Abers2016] or are accessible through web-based applications [@Pommier2011].  In this context, `pide` is a solution fulfilling the need for a versatile library capable of facilitating petrophysical calculations across a range of properties and supporting the creation of specific scientific tools. Beyond this, `pide` aims to host toolkits tailored for specific purposes, such as constructing realistic, petrophysically constrained synthetic models and converting numerical plate tectonic models into synthetic geophysical tomographies.

# Library modules and methods

The general workflow diagram of the library can ben seen in Figure 1. The library has three main classes used in these calculations: `pide`, `material` and `model`. 

`pide` is the main class in which the electrical conductivity and seismic-related observables are calculated. In order to achieve this, the relevant parameters have to be defined in the pide object with the associated functions (e.g., composition, water content, interconnectivity). Just using the `pide` class, for instance, the user can make a figure of all calculations of all the olivine electrical conductivity and seismic velocity models for olivine. While experimental parameters for electrical conductivity parameters are defined and calculated within the `pide`, seismic velocities of given earth materials are calculated through `SAnTex` library automatically.

`material` is the class that can be specified as a holder of pre-defined material properties. For instance, one can create a Lherzolite material by mixing specific modal proportions of ol, opx, cpx and gt, how these constituents are interconnected, or how water behave among them. This class allows the user to define some parameters tied to other parameters via user-defined functions. This allows the users to streamline some specific processes in a modular manner.

The `model` class, on the other hand, is where a collection of `material` objects can be appended with specific positions indexed in 3D space.
`model` is also where the user can calculate the magnetic and gravitational anomalies solely since these observables are dependent on the position of the materials and assigned magnetic and density parameters only. `pide` can generate synthetic data for magnetic and gravitational anomalies utilising the `harmonica` library [@Soler2019]. However, it can only generate synthetic magnetotelluric and seismic velocity models since their forward modelling is not as easy to implement within the pide framework. Currently, pide can generate input files for commonly used magnetotelluric modelling algorithms `ModEM` [@Kelbert2014] and `Mare2DEM` [@Key2016], and then users can generate synthetic data using the algorithms themselves with the given synthetic model.

`pide` also comes with several tools that can exploit the program architecthure. 'model_modifier' functions, for instance, can be used to calculate certain specific hypothesis tests. `inversion` modules, on the other hand, can be utilised to invert for specific input parameters that fit real geophysical data.

![Workflow Chart for pide \label{fig:pide_wflow}](figures/pide_workflow.png)

# Acknowledgements
This study is supported by the Australian Research Council (ARC) Linkage Grant ARC-LP190100146 and ARC DP Grant ARC-DP220100709. 

# References
