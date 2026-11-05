# Migrating from v0.3

Version 0.4 replaces the per-function interface of v0.3 by a single pass over the trajectories
driven by a parameters file. Spectra are no longer read back from `fs_*.h5` files: every
analysis is computed from the spectrum of the current frame, so all analyses of a run use the
same surface settings automatically.

## Script → parameters file

A v0.3/experimental script such as

```julia
fluctuation_spectrum_hmin_allC(; pdb_file, traj_file, output_file=fs_file, lipids,
                                 ref_atoms, L_grid=15.0, atom_mask_radius=5)
lipids_atoms_height(; pdb_file, traj_file, fs_file, output_dir, lipids)
peptide_atoms_height(; pdb_file, traj_files=[traj_file], output_file="PEP_zs.dat", lipids,
                       n_residues=8)
box_dimensions(; traj_file, area_file="A.dat")
area_expansion_modulus(; traj_files=[traj_file], output_file="KA.dat")
peptide_curvature_spectrum(; ..., ref_residue=4, q_max=0.20, q_min=0.01, nqs=19)
peptide_curvature_spectrum_X(; ..., ref_residue=4, ...)
lipids_curvature_spectrum(; ..., q_max=0.20, q_min=0.01, nqs=19)
lipids_thickness_spectrum(; ..., q_max=0.20, q_min=0.01, nqs=19)
peptide_thickness_spectrum(; ..., ref_residue=4, ...)
```

run once per trajectory, followed by block averaging, becomes the parameters file
[`examples/parameters.toml`](https://github.com/amiralih/MembraneAnalysis.jl/blob/main/examples/parameters.toml)
and one command for all trajectories:

```bash
membrane-analysis parameters.toml step5_charmm2gmx.pdb step7_{1..15}.xtc -o results
```

## Function mapping

| v0.3 / experimental | v0.4 |
|:--------------------|:-----|
| `fluctuation_spectrum` | computed in every run; `[surface] fill = "neighbors"`; spectra stored in `timeseries/*.h5` with `save_spectra = true` |
| `fluctuation_spectrum_hmin_allC` | `[surface] mask_radius = …` (all backbone atoms of the proteins in `mask_proteins`) |
| `fluctuation_spectrum_hmin` | `[surface] mask_radius = …`, `mask_residues = [ref_residue]` |
| `box_dimensions`, `area_expansion_modulus` | analysis `box` |
| `fluctuation_hq2_data`, `hq2_analysis` | analysis `undulation` |
| `lipids_atoms_height` | analysis `lipid_heights` |
| `peptide_atoms_height` | analysis `protein_heights` |
| `lipids_sampled_curvature(_X/_Y)` | analysis `lipid_sampled_curvature` (`direction`) |
| `peptide_sampled_curvature` | analysis `protein_sampled_curvature` |
| `lipids_curvature_spectrum` | analysis `lipid_curvature_spectrum` |
| `peptide_curvature_spectrum(_X/_Y)` | analysis `protein_curvature_spectrum` (`direction`) |
| `lipids_thickness_spectrum` | analysis `lipid_thickness_spectrum` |
| `peptide_thickness_spectrum` | analysis `protein_thickness_spectrum` |
| `TCB_analysis(_X/_Y)` | analysis `tcb` (`direction`); plots: see [Analyses](analyses.md#tcb) |
| `find_ref_atoms` | choose `ref_atom` from `<LIPID>_zs.dat` or `<LIPID>_TCB_profile.dat` |
| `curvature`, `thickness` | [`curvature_at`](@ref), [`thickness_at`](@ref) |
| `blocking_error` | [`blocking_error`](@ref) (unchanged) |
| `atom_leaflet`, `atom_leaflet_dynamic` | `MembraneSystem.static_leaflet`, `FrameContext.leaflet` |
| `get_index_pairs` | `MembraneAnalysis.get_index_pairs` |

Removed without replacement: `lipids_density_spectrum`, `peptide_density_spectrum`,
`lipids_radial_distribution`, `peptide_radial_distribution`, `count_interactions` (and the
`tag_files` option of the spectra), `voronoi_shells`, and the project-specific experimental
functions (`retina_*`, `hviz_*`, `membrane_sampled_height`, `syt7_angle`, `height_grid`,
`thickness_grid`, `test_height_dynamic`). They remain available in the v0.3.2 release and the
git history. Custom per-frame analyses can be written with [`foreach_frame`](@ref).

## Behavioural changes

- **Memory:** frames are read into a reused buffer; memory no longer grows with trajectory
  length.
- **Statistics:** results of all trajectories are combined; blocks consist of whole frames
  (v0.3 blocked the list of individual values). Errors of derived quantities (``K_A``,
  ``k_c``) are block-jackknife estimates (v0.3: error propagation; the bootstrap error of
  `hq2_analysis` was always zero because of a bug).
- **Area:** ``\langle L_x L_y\rangle`` instead of ``\langle L_x\rangle\langle L_y\rangle``
  (relative difference ~10⁻⁵).
- **Grid:** one grid for all trajectories of a run, from the mean box of sampled frames.
- **Leaflets:** v0.3 used the leaflets of the previous frame to build the midplane surface;
  v0.4 starts from the leaflets of the structure in every frame and repeats the assignment once
  if any lipid changed leaflet (identical results unless lipids flip-flop).
- **Undulation shells** are limited to modes resolved by the grid (v0.3 could include aliased
  modes on small grids).
- **Protein sign convention:** the leaflet of a protein atom is decided by its height relative
  to the bilayer midplane of the PDB structure (v0.3: ``z > 0``, which assumes a structure
  centred at ``z = 0``), or set explicitly with `leaflet`.
- **Protein residues** are identified by their residue number in the PDB file (v0.3 mixed the
  residue number and the sequential residue index). `protein_sampled_curvature` uses the
  `backbone_atoms` (v0.3: only `CA`).
- **TCB weights** are normalised to sum to one.
- **Plots** are no longer produced by the package (Plots.jl and LaTeXStrings are no longer
  dependencies); profile files contain the data for plotting.
- **Output files** have `#` header lines; read them with `readdlm(file; comments=true)`. The
  per-frame spectra are stored as `nx × ny × n_frames` (v0.3: `n_frames × nx × ny`).
- **Centering threshold:** 20 Å (v0.3.2: 15 Å), configurable with `center_threshold`.
