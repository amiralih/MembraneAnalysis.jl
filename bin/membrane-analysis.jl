#!/usr/bin/env julia
#
# Command-line launcher for MembraneAnalysis.jl (works with every supported Julia version).
#
#   julia --project=/path/to/MembraneAnalysis.jl /path/to/MembraneAnalysis.jl/bin/membrane-analysis.jl \
#         parameters.toml structure.pdb traj1.xtc traj2.xtc ... [-o results] [-n 16]
#
# If MembraneAnalysis is installed in your default environment, `--project=...` can be omitted.
# With Julia ≥ 1.12 you can instead install the `membrane-analysis` app (`pkg> app add
# MembraneAnalysis`) or run `julia -m MembraneAnalysis ...`. Run with --help for all options.

using MembraneAnalysis

exit(MembraneAnalysis.cli_main(ARGS))
