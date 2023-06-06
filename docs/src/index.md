# MembraneAnalysis.jl

`MembraneAnalysis.jl` is a Julia package for analyzing simulations of multi-component lipid bilayers. For membranes simulated in a flat geometry, fluctuation modes of membrane surface height, membrane thickness, and lateral distributions of membrane species can be calculated to be used in determining the mechanical properties of the membrane. The package includes functionality to utilize a novel theoretical framework ([1](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.103.042413), [2](https://pubmed.ncbi.nlm.nih.gov/35927953/)) to determine the spatial extent of the influence of a single lipid molecule (or embedded protein) on the mechanical properties of the surrounding membrane, which allows us to determine curvature and thickness preference of membrane species, as well as spatial correlations between them.

`MembraneAnalysis.jl` is designed to be used by biophysical researchers in need of a fast analysis tool that can be easily built upon to enable implementation of new computational methods to guide theoretical descriptions of membrane physics.

Some of the membrane properties that `MembraneAnalysis.jl` can determine include:

- Bending modulus of the membrane via multiple fluctuation-based approaches
- Area expansion modulus of the membrane
- Relative spontanous curvature of the lipid species
- Nuetral surface of the membrane
- Relative thickness preference of the lipid species

!!! note "Basic Tutorial"

    A tutorial demonstrating an example usage of the package can be found [here](https://github.com/amiralih/MembraneAnalysis.jl/blob/main/tutorial/tutorial.md).

!!! tip "Input Data Organization"

    The files used as input data by the functions in this package include (a subset of):
    
    1. a PDB structure file of the membrane simulation
    2. a trajectory file of the membrane simulation (of any format supported by the `Chemfiles.jl` package), or a list of multiple trajectory files
    3. files containing intermediate values generated as output by other functions in the package (either plain text files with the `.dat` extension, or HDF5 binary data files with `.h5` extension)
    
    We recommend placing the structure file and trajectory file(s) in same directory where the code is run from, and store the generated output files in the same folder as well.

    If you have multiple trajectory files representing consecutive segments of a simulation, make sure that your input list of trajectory files is in the correct order if you wish to use the error estimates for the output values that are calculated using the blocking method ([3](https://link.springer.com/chapter/10.1007/BFb0105461)).

```@contents
Pages = ["fluctuation_analysis.md",
         "curvature_analysis.md",
         "thickness_analysis.md",
         "lipids.md",
         "density_analysis.md",
         "interactions.md",
         "utilities.md"]
```

