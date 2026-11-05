# Getting started

## Installation

```julia
pkg> add MembraneAnalysis
```

To use a development version, clone the repository and work in its directory with
`julia --project=.` (run `pkg> instantiate` once).

## Running an analysis

An analysis needs three kinds of input:

1. a **parameters file** (TOML) describing the lipids, proteins and the analyses to run — start
   from [`examples/parameters.toml`](https://github.com/amiralih/MembraneAnalysis.jl/blob/main/examples/parameters.toml)
   and see [Parameters file](parameters.md);
2. a **PDB file** of the system (e.g. `step5_charmm2gmx.pdb` from CHARMM-GUI); the order of the
   atoms must be the same as in the trajectories;
3. one or more **trajectories** in any format supported by
   [Chemfiles](https://chemfiles.org/chemfiles/latest/formats.html) (`.xtc`, `.trr`, `.dcd`,
   `.nc`, …), given **in chronological order**. Only rectangular boxes are supported.

### From the command line

With Julia ≥ 1.12 the package provides a `membrane-analysis` app
(`pkg> app add MembraneAnalysis`, then make sure `~/.julia/bin` is on your `PATH`); it can also
be called as `julia -m MembraneAnalysis`. With any supported Julia version you can use the
launcher script in the package:

```bash
julia --project=/path/to/MembraneAnalysis.jl /path/to/MembraneAnalysis.jl/bin/membrane-analysis.jl \
      parameters.toml step5_charmm2gmx.pdb step7_{1..15}.xtc -o results
```

Options:

| option | meaning |
|:-------|:--------|
| `-o`, `--output DIR` | output directory (overrides `[output] directory`) |
| `-n`, `--nprocs N` | number of worker processes |
| `--resume` | reuse per-trajectory results of an interrupted run |
| `--chunk-frames N` | split trajectories into work units of at most `N` frames |
| `--check` | validate all inputs, print the plan and exit |
| `--list-analyses` | list the available analysis types |
| `-h`, `--help`, `--version` | help, version |

!!! warning "Order of trajectories"
    Block averaging assumes that the trajectories are given in chronological order. A shell
    glob such as `step7_*.xtc` sorts `step7_10.xtc` before `step7_2.xtc`; use brace
    expansion (`step7_{1..15}.xtc`) or `$(ls -v step7_*.xtc)`. The program warns when the
    order differs from the natural numerical order.

### From Julia

```julia
using MembraneAnalysis
summary = run_analysis("parameters.toml", "step5_charmm2gmx.pdb",
                       ["step7_$(i).xtc" for i in 1:15]; output_dir="results", nprocs=8)
summary["box"]["KA"]
```

## Parallel processing and memory

Each trajectory (or, when there are fewer trajectories than processes, each chunk of a
trajectory) is a *work unit*. Work units are distributed over worker processes started with
Julia's `Distributed` standard library. Every process reads one frame at a time into a reused
buffer, so memory use per process is about the size of one frame plus the (small) per-frame
results.

The number of processes is, in order of precedence: `--nprocs`, the CPUs allocated by the batch
system (`SLURM_CPUS_PER_TASK`, `SLURM_CPUS_ON_NODE`, `PBS_NP`, `NSLOTS`), or all CPU threads of
the machine. If Julia was started with workers already (`julia -p N` or
`julia --machine-file hosts`), those are used instead; this is how several nodes can be used.

On cluster file systems without file locking (some Lustre/NFS setups) HDF5 may fail to create
files; set `export HDF5_USE_FILE_LOCKING=FALSE` in the job script in that case.

An example SLURM script is in
[`examples/run_analysis.slurm`](https://github.com/amiralih/MembraneAnalysis.jl/blob/main/examples/run_analysis.slurm).
Use `--resume` in job scripts: if a job is killed, rerunning it only processes the trajectories
that were not finished.

## Output

```
results/
├── summary.toml              scalar results (A, K_A, k_c, neutral surface, …) and run info
├── parameters_used.toml      the parameters with the defaults of the global sections filled in
├── A.dat, KA.dat, …          one file per result (see "Analyses and output files")
└── timeseries/
    ├── 001_step7_1.h5        per-frame values of every analysis for trajectory 1
    └── …
```

Result files are tab-separated text with `#` header lines (run information, description, column
names). Read them with

```julia
using DelimitedFiles
data = readdlm("results/POPC_cqs.dat"; comments=true)
```

The HDF5 files contain, for every analysis, a group `analyses/<id>` with datasets `values`
(per-frame value of every component, `n_components × n_frames`), `counts` (number of samples
behind every value), `sums` (`values .* counts`, used for exact resuming) and the attribute
`labels`; `box` (`3 × n_frames`); and, if
`save_spectra = true`, `spectra/hq` and `spectra/tq` (`nx × ny × n_frames`, complex). They can
be used to recompute statistics in any way (from Python, `h5py` shows the dimensions in reverse
order).
