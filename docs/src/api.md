# Julia API

## Running analyses

```@docs
run_analysis
load_parameters
parse_parameters
AnalysisConfig
ParameterError
default_nprocs
MembraneAnalysis.cli_main
```

## System definition

```@docs
Lipid
Protein
MembraneSystem
read_pdb
PDBAtoms
```

Lipid presets: `POPC_aa`, `Chol_aa`, `DMPC_aa`, `DOPC_aa`, `SOPC_aa` (CHARMM36 all-atom) and
`POPC_m2`, `DOPE_m2`, `DOPC_m2` (Martini 2).

## Building blocks for custom analyses

```@docs
foreach_frame
FrameContext
TrajectoryReader
with_trajectory
read_frame!
prepare_frame!
SurfaceSettings
SurfaceWorkspace
compute_surface!
grid_size_for
QBins
spectrum_bins
q_shell
FourierField
evaluate!
curvature_at
thickness_at
```

## Statistics

```@docs
blocking_error
equal_blocks
trajectory_blocks
optimal_block_length
BlockStats
block_means
block_sems
jackknife
```

## Analyses

Analyses are normally created from a parameters file; the types can also be constructed
directly and are documented here for reference.

```@docs
AbstractAnalysis
BoxAnalysis
UndulationAnalysis
LipidHeights
ProteinHeights
LipidSampledCurvature
ProteinSampledCurvature
LipidSpectrum
ProteinSpectrum
TCBAnalysis
```

### Writing a new analysis

An analysis is a struct `<: AbstractAnalysis` with methods for

- `MembraneAnalysis.analysis_id(a)` — a unique string;
- `MembraneAnalysis.component_labels(a)` — labels of the per-frame components;
- `MembraneAnalysis.output_files(a, system)` — names of the files that `write_results` writes;
- `MembraneAnalysis.accumulate_frame!(sums, counts, a, system, ctx)` — add the contributions of
  one frame (`ctx::FrameContext`) to the component vectors `sums` and `counts`;
- `MembraneAnalysis.write_results(a, stats::BlockStats, system, out)` — write result files
  (with `MembraneAnalysis.write_table`) and return a dictionary of scalar results;
- optionally `MembraneAnalysis.needs_surface(a)` (`true` by default).

To make it available in parameters files, add its keys to `_ANALYSIS_KEYS` and a branch to
`_build_analysis` in `src/parameters.jl`.
