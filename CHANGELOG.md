# Changelog

All notable changes to this project are documented in this file. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[semantic versioning](https://semver.org/) (before 1.0, minor versions may break the API).

## [0.4.0] — 2026-11-05

This release is a redesign of the package. See the
[migration guide](https://amiralih.github.io/MembraneAnalysis.jl/dev/migration/).

### Added

- Command-line program `membrane-analysis` (Julia ≥ 1.12 app, `julia -m MembraneAnalysis`, or
  `bin/membrane-analysis.jl`) and [`run_analysis`](https://amiralih.github.io/MembraneAnalysis.jl/dev/api/):
  runs the analyses defined in a TOML parameters file on all trajectories.
- Parallel processing of trajectories (and chunks of trajectories) with `Distributed`; the
  number of processes defaults to the CPUs allocated by SLURM/PBS/SGE.
- Final statistics over all trajectories: means, block-averaged standard errors (`equal`,
  `trajectory` or automatic block size) and block-jackknife errors of derived quantities.
- Per-trajectory HDF5 time series of every analysis; `--resume` reuses them after an
  interruption.
- Reconstruction of the membrane surface around proteins by bending-energy minimisation
  (`[surface] mask_radius`, `fill = "min_curvature"`), using a kernel formulation that needs one
  FFT pair per reconstruction.
- Directional (`direction = "x"`/`"y"`) curvature spectra, sampled curvature and TCB.
- Thickness spectra of proteins, thickness undulation spectrum (`tq2.dat`).
- `Protein` definitions (segments, chains or segment prefix; backbone atoms; leaflet sign).
- Validation of lipid/protein definitions against the structure and of parameters files.
- `foreach_frame` and other building blocks for custom analyses.

### Changed

- Trajectories are read into a single reused `Chemfiles.Frame`: memory use no longer grows with
  the number of frames.
- Each frame is read once and all analyses are computed in the same pass; spectra are not
  re-read from `fs_*.h5` files and the PDB file is parsed once (previously once per frame for
  the dynamic leaflet assignment).
- Blocks for error estimates consist of whole frames; the area is ⟨Lx Ly⟩; TCB weights are
  normalised; protein leaflet signs are determined relative to the bilayer midplane of the
  structure; protein residues are identified by PDB residue number.
- The leaflet assignment of a frame no longer depends on the previous frame: it starts from
  the leaflets of the structure and is refined once if lipids changed leaflet, so results do not
  depend on how trajectories are split between processes.
- Output text files have `#` header lines. Spectra in HDF5 files are stored as
  `nx × ny × n_frames`.
- Minimum Julia version 1.10. Dependencies PDBTools, Plots, LaTeXStrings, Distances, StatsBase
  and ProgressMeter were removed (a small built-in PDB reader is used).

### Fixed

- The parametric-bootstrap error of the undulation bending modulus was always zero (the fits
  used the unperturbed data).
- The `HC` undulation model returned a scalar instead of a vector.
- Atoms exactly at the upper box edge could produce an out-of-range grid cell.
- Undulation shells could include aliased modes (indices beyond the Nyquist frequency) on small
  grids; the number of shells is now limited to the resolved modes.
- Residue numbers with five digits (CHARMM-GUI PDB files) are read correctly.
- The z-centering loop could run forever for non-flat membranes; it now stops with an error.

### Removed

- The v0.3 functions (`fluctuation_spectrum`, `lipids_curvature_spectrum`, …), replaced by
  analyses in the parameters file.
- Density spectra, radial distribution functions, interaction counting and tag filtering.
- Plotting (`tcb_plot`, `hq2_plot`); profile files contain the data instead.
- Experimental, project-specific functions of the development branch.

## [0.3.2]

Last release of the v0.3 series (see the git history).
