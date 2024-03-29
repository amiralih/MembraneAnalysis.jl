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
 - name: Eunice Kennedy Shriver National Institute of Child Health and Human Development, Bethesda, MD, United States of America
   index: 1
date: 1 December 2022
bibliography: paper.bib

---

# Summary

Biological membranes separate living cells from their surroundings, and, in case of eukaryotes, partition the cell into its constituent organelles. These bilayer structures, which consist of hundreds of different lipid species and associated proteins, play a central role in a variety of vital biological processes. Molecular dynamics (MD) simulations are a potent tool, used in conjunction with theoretical modeling and experimental studies, to investigate the physical and chemical properties of biomembranes. A crucial step in leveraging the power of MD in biophysics is analyzing the data produced by simulations in order to extract equilibrium and dynamic quantities of interest, which will inform and validate our theoretical models and help interpret experimental results. Availability of computationally efficient, flexible, and extensible software facilitates researchers in this endeavor.

# Statement of need

`MembraneAnalysis.jl` is a Julia package for analyzing simulations of multi-component lipid bilayers. For membranes simulated in a flat geometry, fluctuation modes of membrane surface height, membrane thickness, and lateral distributions of membrane species can be calculated to be used in determining the mechanical properties of the membrane. The package includes functionality to utilize a novel theoretical framework [@sapp2021spatial; @lessen2022molecular] to determine the spatial extent of the influence of a single lipid molecule (or embedded protein) on the mechanical properties of the surrounding membrane, which allows us to determine curvature and thickness preference of membrane species, as well as spatial correlations between them.

`MembraneAnalysis.jl` is designed to be used by biophysical researchers in need of a fast analysis tool that can be easily built upon to enable implementation of new computational methods to guide theoretical descriptions of membrane physics.

Some of the membrane properties that `MembraneAnalysis.jl` can determine include:

- Bending modulus of the membrane via multiple fluctuation-based approaches
- Area expansion modulus of the membrane
- Relative spontaneous curvature of the lipid species
- Neutral surface of the membrane
- Relative thickness preference of the lipid species

# State of the field

There are quite a few software programs and packages developed for processing molecular dynamics simulations of lipid membranes, which can perform a number of analyses. Tools such as GridMAT-MD [@allen2009gridmat], APL@ Voro [@lukat2013apl], MEMBPlugin [@guixa2014membplugin], FATSLiM [@buchoux2017fatslim], MemSurfer [@bhatia2019memsurfer], and LiPyphilic [@smith2021lipyphilic] calculated properties like membrane thickness, area-per-lipid, and order parameter. [MembraneCurvature](https://github.com/MDAnalysis/membrane-curvature) is an [MDAnalysis](https://www.mdanalysis.org/) tool to calculate membrane curvature. What sets `MembraneAnalysis.jl` apart from these tools is providing functionality to calculate various additional physical properties of interest, such as those  describing intrinsic curvature and thickness preferences of the different lipid types present in the membrane, and the moduli characterizing the elastic properties of the membrane. These quantities are crucial in understanding any biological process that involves reshaping of the membrane.

# Acknowledgements

This project was supported by the Intramural Research Program (IRP) of the Eunice Kennedy Shriver National Institute of Child Health and Human Development (NICHD) at the National Institutes of Health (NIH).

# References
