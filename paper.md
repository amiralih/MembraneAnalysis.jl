---
title: 'MembraneAnalysis.jl: A Julia package for analyzing molecular dynamics simulations of lipid membranes'
tags:
  - Julia
  - biophysics
  - molecular dynamics
  - lipid membranes
authors:
  - name: Amirali Hossein
    orcid: 0000-0002-2580-3577
    affiliation: 1
  - name: Alexander J. Sodt
    orcid: 0000-0002-5570-8212
    affiliation: 1
affiliations:
 - name: Eunice Kennedy Shriver National Institute of Child Health and Human Development, Bethesda, MD 20892, United States
   index: 1
date: 1 December 2022
bibliography: paper.bib

---

# Summary

Biological membranes separate living cells from their surrounding, and, in case of eukaryotes, partition the cell into its constituent organelles. These bilayer structures, which consist of hundreds of different lipid species and associated proteins, play a central role in a variety of vital biological processes. Molecular dynamics (MD) simulations are a potent tool, used in conjunction with theoretical modeling and experimental studies, to investigate the physical and chemical properties of biomembranes. A crucial step in leveraging the power of MD in biophysics is analyzing the data produced by simulations in order to extract equilibrium and dynamic quantities of interest, which will inform and validate our theoretical models and help interpret experimental results. Availability of computationally efficient, flexible, and extensible software facilitates researchers in this endeavor.

# Statement of need

`MembraneAnalysis.jl` is a Julia package for analyzing simulations of multi-component lipid bilayers. For membranes simulated in a flat geometry, fluctuation modes of membrane surface height, membrane thickness, and lateral distributions of membrane species can be calculated to be used in determining the mechanical properties of the membrane. The package includes functionality to utilize a novel theoretical framework [@sapp2021spatial; @lessen2022molecular] to determine the spatial extent of the influence of a single lipid molecule (or embedded protein) on the mechanical properties of the surrounding membrane, which allows us to determine curvature and thickness preference of membrane species, as well as spatial correlations between them.

`MembraneAnalysis.jl` is designed to be used by biophysical researchers in need of a fast analysis tool that can be easily built upon to enable implementation of new computational methods to guide theoretical descriptions of membrane physics.

# Acknowledgements

This project was supported by the Intramural Research Program (IRP) of the Eunice Kennedy Shriver National Institute of Child Health and Human Development (NICHD) at the National Institutes of Health (NIH).

# References
