# Parameters file

The parameters file is a [TOML](https://toml.io) file. Unknown keys are reported as errors, so
typos do not go unnoticed. A complete annotated example is
[`examples/parameters.toml`](https://github.com/amiralih/MembraneAnalysis.jl/blob/main/examples/parameters.toml).
`membrane-analysis --check …` validates a parameters file against a structure and trajectories
without running the analysis, and every run writes `parameters_used.toml` with the defaults of
the global sections filled in (the resolved grid size is recorded in `summary.toml`).

## Top level

| key | default | description |
|:----|:--------|:------------|
| `title` | `""` | free text shown in the log |

## `[[lipids]]` (required, one table per species)

| key | default | description |
|:----|:--------|:------------|
| `name` | required | residue name in the PDB file |
| `head_atom` | required | first atom of every molecule; its height relative to the midplane surface decides the leaflet |
| `ref_atom` | required | reference atom: defines the membrane surface and the lipid position in spectra |
| `tail_atom` | required | tail-end atom; the mean height of all tail atoms is the bilayer midplane |
| `n_atoms` | required | number of atoms per molecule (atoms must be contiguous, starting with `head_atom`) |
| `height_atoms` | atoms named `C*`, `D*`, `G*`, `X*` | atoms used by `lipid_heights`, `lipid_sampled_curvature` and `tcb` |
| `preset` | — | start from a predefined lipid (`POPC_aa`, `Chol_aa`, `DMPC_aa`, `DOPC_aa`, `SOPC_aa`, `POPC_m2`, `DOPE_m2`, `DOPC_m2`); other keys override its fields |

The definitions are checked against the structure: missing atoms or an `n_atoms` that does not
match the residues raise an error.

## `[[proteins]]` (optional)

| key | default | description |
|:----|:--------|:------------|
| `name` | required | label used in output file names |
| `segments` | — | exact segment names (PDB columns 73–76) |
| `chains` | — | chain identifiers (used if `segments` is not given) |
| `segment_prefix` | `"PRO"` | segment-name prefix (used if neither `segments` nor `chains` is given) |
| `backbone_atoms` | `["BB", "CA"]` | atom names representing a residue |
| `leaflet` | `"auto"` | sign convention for curvature and heights: `"upper"` (+1), `"lower"` (−1), or `"auto"`: +1 for atoms above the bilayer midplane in the PDB structure, −1 below |

Residues are identified by their residue number in the PDB file. All copies selected by one
protein definition (e.g. segments `PROA` and `PROB`) are averaged in per-residue results.

## `[surface]`

| key | default | description |
|:----|:--------|:------------|
| `grid_spacing` | `15.0` | target cell size in Å; the grid has `ceil(⟨L⟩ / grid_spacing)` cells per direction, with `⟨L⟩` estimated from evenly spaced frames of every trajectory (up to 50 per trajectory, about 500 in total) |
| `grid_size` | — | number of cells `[nx, ny]` (overrides `grid_spacing`) |
| `mask_radius` | `0.0` | cells whose centre is closer than this (in the xy plane, Å) to a backbone atom of a masked protein are reconstructed by bending-energy minimisation |
| `mask_proteins` | all proteins if `mask_radius > 0` | proteins whose backbone atoms define the mask |
| `mask_residues` | all | restrict the mask to these residues (list or range such as `"1:8"`) |
| `fill` | `"min_curvature"` with a mask, otherwise `"neighbors"` | treatment of cells without reference atoms: `"neighbors"` (sum of the four neighbouring cells, v0.3 behaviour) or `"min_curvature"` (bending-energy minimisation); a mask requires `"min_curvature"` |
| `center_threshold` | `20.0` | Å; if the standard deviation of tail-atom heights exceeds this, the bilayer is assumed to be split across the periodic boundary in z and the system is shifted |

See [Methods](methods.md) for details.

## `[spectrum]`

Defaults for all `*_spectrum` analyses; every analysis may override them.

| key | default | description |
|:----|:--------|:------------|
| `q_min` | `0.01` | lower end of the first bin (Å⁻¹) |
| `q_max` | `0.21` | upper end of the last bin (Å⁻¹) |
| `n_bins` | `10` | number of equal-width bins |

## `[statistics]`

| key | default | description |
|:----|:--------|:------------|
| `blocking` | `"equal"` | `"equal"`: `n_blocks` contiguous blocks of `⌊N/n_blocks⌋` frames of the concatenated trajectories; `"trajectory"`: one block per trajectory file; `"auto"`: block length chosen automatically for every analysis |
| `n_blocks` | `10` | number of blocks for `"equal"` |

## `[output]`

| key | default | description |
|:----|:--------|:------------|
| `directory` | `"membrane_analysis_results"` | output directory (the `-o` option takes precedence) |
| `save_spectra` | `true` | store the height and thickness spectra of every frame in the time-series files |

## `[[analysis]]` (required, one table per analysis)

Every analysis table has a `type` and an optional `label`, which replaces the default file-name
suffix (files are named `<LIPID>_<label>.dat` or `<PROTEIN>_<label>.dat`; for `box`,
`undulation` and `tcb` the label is appended to the file names). Two analyses of the same type
with the same label (for the same lipids), or two analyses that would write the same file, are
rejected.

| `type` | keys | output |
|:-------|:-----|:-------|
| `box` | — | `box_dims.dat`, `A.dat`, `KA.dat` |
| `undulation` | `q_max` (required), `model` (`"TD"`), `n_fit` (all) | `hq2.dat`, `tq2.dat`, `kc_hq2.dat` |
| `lipid_heights` | `lipids` (all) | `<LIPID>_zs.dat` |
| `lipid_sampled_curvature` | `q_max` (required), `q_min` (0), `direction`, `lipids` | `<LIPID>_cs.dat` |
| `lipid_curvature_spectrum` | `q_min`, `q_max`, `n_bins`, `direction`, `lipids`, `ref_atoms` | `<LIPID>_cqs.dat` |
| `lipid_thickness_spectrum` | same as above | `<LIPID>_tqs.dat` |
| `protein_heights` | `protein`, `residues` (all) | `<PROTEIN>_zs.dat` |
| `protein_sampled_curvature` | `protein`, `q_max` (required), `q_min` (0), `direction`, `residues` | `<PROTEIN>_cs.dat` |
| `protein_curvature_spectrum` | `protein`, `residue` (required), `q_min`, `q_max`, `n_bins`, `direction` | `<PROTEIN>_cqs.dat` |
| `protein_thickness_spectrum` | same as above | `<PROTEIN>_tqs.dat` |
| `tcb` | `q_max` (required), `z_cutoff` (required), `q_min` (0), `weights` (equal), `direction`, `lipids` | `TCB.dat`, `TCB_lipids.dat`, `<LIPID>_TCB_profile.dat` |

- `direction`: `"all"` (default), `"x"` or `"y"` — use only Fourier modes with wave vectors
  along x or y. Directional analyses get `_X`/`_Y` appended to their default file names.
- `lipids`: list of lipid names to analyse (default: all).
- `protein`: protein name; may be omitted when only one protein is defined.
- `ref_atoms`: table of reference atoms replacing the lipids' `ref_atom`, e.g.
  `ref_atoms = { POPC = "D2A" }`.
- `weights`: area fractions of the lipids for TCB, e.g. `weights = { DOPC = 0.9, CHL1 = 0.1 }`
  (normalised to sum 1); if given, every analysed lipid needs a weight.
- `model`: `"HC"` (Helfrich–Canham), `"MNK"` (May–Narang–Kopelevich) or `"TD"`
  (Terzi–Deserno).

See [Analyses and output files](analyses.md) for definitions.
