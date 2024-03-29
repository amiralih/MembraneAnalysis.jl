# Example usage of `MembraneAnalysis.jl`

In this tutorial we showcase some of the main functionality of the package by analyzing an atomistic simulations of a membrane containing DOPC lipids and %10 mol fraction cholesterol. `MembraneAnalysis.jl` needs a PDB structure file (here `s.pdb`) and trajectory files of any format supported by `Chemfiles.jl` (here `t1.dcd` to `t10.dcd`). The files needed for this tutorial can be found [here](https://doi.org/10.5281/zenodo.8007968). To reproduce the analysis in this tutorial, first unzip the compressed folder and navigate into it, and then run the following Julia code snippets in order.

## Setting analysis parameters

It is useful to define a set of parameters that we will need for input arguments of the methods.
- `pdb_file`: the PDB structure file
- `traj_dir`: folder containing trajectory files
- `traj_name`: trajecotry file name without numerical index and file extension
- `traj_inds`: range of numerical indices of trajectory files
- `traj_ext`: file extension of trajecotries
- `lipids`: list of lipids in the system (from instances of `Lipid` type defined in `src/lipids.jl`
- `L_grid`: lattice grid length (in length unit of trajectories) used for discretizing the membrane surface
- `q_max`: q mode magnitude upper cut-off in the analysis
- `output_dir`: a directory for all the output files generated by methods

```julia
using MembraneAnalysis

begin
  const pdb_file = "s.pdb"
  const traj_dir = "./"
  const traj_name = "t"
  const traj_inds = 1:10
  const traj_ext = ".dcd"
  const lipids = [DOPC_aa, Chol_aa]
  const L_grid = 15
  const q_max = 0.08
  const output_dir = "./"
  if !isdir(output_dir) mkdir(output_dir) end
end
```

## Calculating fluctuation spectrum

`fluctuation_spectrum` calculates discrete Fourirer modes of surface height and thickness fluctuations and stores them in an HDF5 file for each trajectory.

```julia
for i in traj_inds
    traj_file = traj_dir * traj_name * string(i) * traj_ext
    output_file = output_dir * "fs_$(i).h5"

    fluctuation_spectrum(;
        pdb_file=pdb_file,
        traj_file=traj_file,
        output_file=output_file,
        lipids=lipids,
        L_grid=L_grid
    )
end
```

## Calculating area expansion modulus

`box_dimensions` calculates the area and stores it in a file.

```julia
traj_file = traj_dir * traj_name * string(first(traj_inds)) * traj_ext
output_file = output_dir * "KA.dat"

box_dimensions(;
    traj_file=traj_file,
    area_file=output_dir * "A.dat"
)
```

`area_expansion_modulus` calculates area expansion modulus from area fluctuations and stores them in a file.

```julia
traj_files = [traj_dir * traj_name * string(i) * traj_ext for i in traj_inds]
output_file = output_dir * "KA.dat"

area_expansion_modulus(;
    traj_files=traj_files,
    output_file=output_file
)
```

We can find the calculated value of $\approx 0.7 \frac{k_{\mathrm{B}}T}{\mathrm{Å}^2}$ for area expansion modulus of the membrane in `KA.dat`.

## Calculating mean height from midplane of the heavy atoms of the lipids

We use the first trajecotry to calculate the average height of the heavy atoms of the lipids using `lipids_atoms_heights` method which will save them to `XXXX_zs.dat` files in the specified output directory (where "XXXX" is the name of the lipid).

```julia
traj_file = traj_dir * traj_name * string(first(traj_inds)) * traj_ext
fs_file = output_dir * "fs_$(first(traj_inds)).h5"

lipids_atoms_height(;
    pdb_file=pdb_file,
    traj_file=traj_file,
    fs_file=fs_file,
    output_dir=output_dir,
    lipids=lipids
)
```

## Calculating mean sampled curvature of heavy atoms of the lipids

We calculate the average sampled curvature of the heavy atoms of the lipids using `lipids_sampled_curvature` method which will save them to `XXXX_cs.dat` files in the specified output directory (where "XXXX" is the name of the lipid).

```julia
traj_files = [traj_dir * traj_name * string(i) * traj_ext for i in traj_inds]
fs_files = [output_dir * "fs_$(i).h5" for i in traj_inds]

lipids_sampled_curvature(;
    pdb_file=pdb_file,
    traj_files=traj_files,
    fs_files=fs_files,
    output_dir=output_dir,
    lipids=lipids,
    q_max=q_max
)
```

## Transverse curvature bias analysis

Using the outputs of the two previous steps we can plot the transverse curvature bias curves, which can be used to identify the neutral surface atom of each lipid, calculate the bending modulus of the membrane and estimate the difference in spontanous curvature of the membrane components [[1]](#1).

```julia
TCB_analysis(;
    input_dir=output_dir,
    lipids=lipids,
    weights=[0.9, 0.1],
    z_cutoff=14,
    output_dir=output_dir,
    tcb_plot=true,
)
```

<p align="center">
<img src="https://github.com/amiralih/MembraneAnalysis.jl/blob/13e6a3dc86be379ad39753cf3d9dfe4effc3c919/tutorial/TCB_DOPC_10.png" width="500">
</p>

We can find the calculated value of $\approx 20 k_{\mathrm{B}}T$ for bending modulus of the membrane in `kc.dat`.

## References
<a id="1">[1]</a> 
Sapp, K. C., Beaven, A. H., and Sodt, A. J. (2021). 
Spatial extent of a single lipid's influence on bilayer mechanics
Phys. Rev. E 103, 042413.
