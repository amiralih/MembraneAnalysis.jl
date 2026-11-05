# MembraneAnalysis.jl

```@raw html
<p align="center"><img src="assets/logo.png" width="400"></p>
```

`MembraneAnalysis.jl` is a Julia package for analysing molecular dynamics simulations of lipid
bilayers simulated in a flat geometry. From the trajectories it computes the Fourier modes of
membrane height and thickness fluctuations and uses them to determine mechanical properties of
the membrane and the curvature and thickness preferences of lipids and embedded proteins
([Sapp, Beaven & Sodt, *Phys. Rev. E* 103, 042413 (2021)](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.103.042413);
[Hossein & Sodt, *J. Open Source Softw.* 8, 5380 (2023)](https://doi.org/10.21105/joss.05380)).

Properties that can be determined include:

- the bending modulus of the membrane, from transverse curvature bias (TCB) and from the
  undulation spectrum,
- the area expansion modulus,
- curvature and thickness spectra of lipids and of proteins/peptides (optionally restricted
  to modes along x or y),
- the relative spontaneous curvature of lipid species and the neutral surface,
- heights of lipid atoms and protein residues relative to the bilayer midplane.

Around embedded proteins the membrane surface can be reconstructed by minimising its bending
energy, so that the protein itself does not distort the height and thickness spectra.

## How it works

You describe the system (lipids, proteins) and the analyses you want in a
[parameters file](parameters.md) and run

```bash
membrane-analysis parameters.toml structure.pdb traj_1.xtc traj_2.xtc ... -o results
```

All trajectories are processed in parallel, using as many processes as CPUs are available;
each process holds a single frame in memory at a time. The per-frame results of every
trajectory are saved, and the final means and standard errors are obtained by block averaging
over the whole (concatenated) simulation. See [Getting started](getting_started.md).

!!! note "Version 0.4"
    Version 0.4 is a redesign of the package with a new, configuration-file based interface.
    Functions of v0.3 such as `fluctuation_spectrum` or `lipids_curvature_spectrum` have been
    replaced; see [Migrating from v0.3](migration.md).

```@contents
Pages = ["getting_started.md", "parameters.md", "analyses.md", "methods.md", "api.md", "migration.md"]
Depth = 1
```
