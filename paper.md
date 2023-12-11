---
#title: 'pide: an Earth material library for petrophysical calculations and synthetic modelling of Earth systems.'
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
  - name: Utpal Singh
    affiliation: 1
  - name: Lu Li
    affiliation: 2
  - name: Patrice F. Rey
    affiliation: 1
  - name: Maria Constanza Manassero
    affiliation: 3

affiliations:
 - name: School of Geosciences, University of Sydney, Sydney, Australia.
   index: 1
 - name: School of Earth Sciences, University of Western Australia, Perth, Australia.
   index: 2
 - name: School of Natural Sciences, University of Tasmania, Hobart, Australia.
   index: 3
date: 16 August 2023
bibliography: paper.bib
---

# Summary

# Statement of need
Given the inherent heterogeneity and complexity of Earth systems, geophysical tomographies often yield complex 2D and 3D images that are challenging to interpret. To enhance their interpretations, researchers commonly turn to experimental petrology and mineral physics, covering various geophysical properties, including electrical conductivity, magnetic susceptibility, seismic velocity, and rheology. These properties are sensitive to phase transitions, partial melting, major and trace elements partitioning, mineral solubilities, and phase-mixing models. Numerous specialized tools have been designed to address specific properties, many of which feature graphical interfaces [@Ozaydin2020; @Abers2016] or are accessible through web-based applications [@Pommier2011].  In this context, `pide` is a solution fulfilling the need for a versatile library capable of facilitating petrophysical calculations across a range of properties and supporting the creation of specific scientific tools. Beyond this, `pide` aims to host toolkits tailored for specific purposes, such as constructing realistic, petrophysically constrained synthetic models and converting numerical plate tectonic models into synthetic geophysical tomographies.

# Package Summary

`pide` is a Python3 library to calculate geophysical observables employing the results from experimental petrology, mineral/rock physics and thermodynamic modelling studies. The geophysical observables currently included in the library are electrical conductivity, seismic velocity, seismic anisotropy, and magnetic and gravitational anomalies.

The general workflow diagram of the library can ben seen in Figure 1. The library has three main classes used in these calculations: `pide`, `material` and `model`. 

`pide` is the main class in which the electrical conductivity and seismic-related observables are calculated. In order to achieve this, the relevant parameters have to be defined in the SEL object with the associated functions (e.g., composition, water content, interconnectivity). Just using the `pide` class, for instance, the user can make a figure of all calculations of all the olivine electrical conductivity and seismic velocity models for olivine (Figure 2). While experimental parameters for electrical conductivity and magnetic properties are defined and calculated within the `pide`, seismic velocity and anisotropy can be calculated by @pyUtpal library.

`material` is the class that can be specified as a holder of pre-defined material properties. For instance, one can create a Lherzolite material by mixing specific modal proportions of ol, opx, cpx and gt. This class allows the user to define some parameters tied to other parameters via user-defined functions. For example, water content can be defined as a function with changing modal opx and temperature in the material.

`model` class, on the other hand, is where a collection of materials can be appended with specific positions in 3D space. This `model` field can be associated with the external mapping of `material` objects. Figure 3, we can see the calculated geophysical observables (electrical conductivity and seismic velocity) of a model where the material field is extracted from a numerical model of a xxx zone.

The `model` class, on the other hand, is where a collection of `materials`  can be appended with specific positions in 3D space. This `model` field can be associated with the external mapping of `material` objects. In Figure 3, we can see the calculated geophysical observables (electrical conductivity and seismic velocity) of a `model` where the material field is extracted from a numerical model of a xxx zone. `model` is also where the user can calculate the magnetic and gravitational anomalies solely since these observables are dependent on the position of the materials and instantly assigned magnetic and density parameters only. 

`pide` also comes with several tools that can exploit this structure. 'model_modifier' functions, for instance, can be used to calculate certain specific hypothesis tests. `inversion` modules, on the other hand, can be utilised to invert for specific input parameters that fit real geophysical data.

![Workflow Chart for pide \label{fig:pide_wflow}](SEL/docs/figures/pide_workflow.png)

# Acknowledgements
This study is supported by the Australian Research Council (ARC)Linkage Grant #Grantnumber and ARC DP Grant #Grantnumber. 

# References
