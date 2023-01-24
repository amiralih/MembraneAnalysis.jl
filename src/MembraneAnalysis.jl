module MembraneAnalysis

import PDBTools
import Chemfiles
using Statistics
using FFTW
using ProgressMeter
using HDF5
using DelimitedFiles
using StatsBase
using Distances
using Plots
using LsqFit

export Lipid,
       POPC_aa,
       Chol_aa,
       DMPC_aa,
       DOPC_aa,
       SOPC_aa,
       POPC_m2,
       DOPE_m2,
       DOPC_m2
include("lipids.jl")

export atom_leaflet,
       atom_leaflet_dynamic,
       box_dimensions,
       blocking_error,
       lipids_atoms_height,
       peptide_atoms_height,
       get_index_pairs,
       find_ref_atoms
include("utilities.jl")

export fluctuation_spectrum,
       area_expansion_modulus
include("fluctuation_analysis.jl")

export curvature,
       lipids_sampled_curvature,
       peptide_sampled_curvature,
       lipids_curvature_spectrum,
       peptide_curvature_spectrum
include("curvature_analysis.jl")

export lipids_thickness_spectrum,
       peptide_thickness_spectrum
include("thickness_analysis.jl")

export lipids_density_spectrum,
       peptide_density_spectrum,
       lipids_radial_distribution,
       peptide_radial_distribution
include("density_analysis.jl")

end # module
