"""
    MembraneAnalysis

Analysis of molecular dynamics simulations of lipid membranes: height and thickness
fluctuation spectra, curvature and thickness spectra of lipids and proteins, transverse
curvature bias, bending and area expansion moduli.

The main entry point is [`run_analysis`](@ref) (or the `membrane-analysis` command), which
runs the analyses defined in a TOML parameters file on a set of trajectories in parallel and
writes means and block-averaged standard errors. See the documentation for details.
"""
module MembraneAnalysis

using Distributed
using FFTW
using LinearAlgebra
using Printf
using Statistics
import Chemfiles
import Dates
import HDF5
import LsqFit
import SHA
import TOML

include_dependency(joinpath(dirname(@__DIR__), "Project.toml"))

"Version of the MembraneAnalysis package."
const PACKAGE_VERSION = something(pkgversion(@__MODULE__),
    VersionNumber(TOML.parsefile(joinpath(dirname(@__DIR__), "Project.toml"))["version"]))

include("lipids.jl")
include("proteins.jl")
include("pdb.jl")
include("system.jl")
include("trajectory.jl")
include("surface.jl")
include("fourier.jl")
include("statistics.jl")
include("output.jl")
include("analyses/interface.jl")
include("analyses/membrane.jl")
include("analyses/local.jl")
include("analyses/tcb.jl")
include("parameters.jl")
include("pipeline.jl")
include("cli.jl")

# running analyses
export run_analysis, load_parameters, parse_parameters, AnalysisConfig, ParameterError,
       AnalysisError, default_nprocs

# system definition
export Lipid, Protein, MembraneSystem, read_pdb, PDBAtoms,
       POPC_aa, Chol_aa, DMPC_aa, DOPC_aa, SOPC_aa, POPC_m2, DOPE_m2, DOPC_m2

# building blocks for custom analyses
export TrajectoryReader, with_trajectory, read_frame!, prepare_frame!, foreach_frame,
       FrameContext, SurfaceSettings, SurfaceWorkspace, compute_surface!, grid_size_for,
       QBins, spectrum_bins, q_shell, FourierField, evaluate!, curvature_at, thickness_at

# statistics
export blocking_error, equal_blocks, trajectory_blocks, optimal_block_length, BlockStats,
       block_means, block_sems, jackknife

# analyses
export AbstractAnalysis, BoxAnalysis, UndulationAnalysis, LipidHeights, ProteinHeights,
       LipidSampledCurvature, ProteinSampledCurvature, LipidSpectrum, ProteinSpectrum,
       TCBAnalysis

end # module
