#!/bin/bash
# Example SLURM job: analyse all trajectories of one replica on one node.
#
# Submit with:  sbatch run_analysis.slurm
# The number of worker processes is taken from SLURM_CPUS_PER_TASK automatically.

#SBATCH --job-name=membrane-analysis
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00

module load julia          # adapt to your cluster

PKG=/data/$USER/MembraneAnalysis.jl          # package directory (or omit --project if installed)
SYSTEM=/data/$USER/aah_anisotropy_project/l8_rep1

cd "$SYSTEM"
# export HDF5_USE_FILE_LOCKING=FALSE   # needed on some Lustre/NFS file systems

# Trajectories must be listed in chronological order: brace expansion keeps the numeric order
# (a glob such as step7_*.xtc would put step7_10 before step7_2).
julia --project="$PKG" "$PKG/bin/membrane-analysis.jl" \
    parameters.toml step5_charmm2gmx.pdb step7_{1..15}.xtc \
    --output results_R5 \
    --resume                 # if the job is restarted, finished trajectories are reused

# Multi-node runs: start Julia with a machine file (or use SlurmClusterManager.jl), e.g.
#   julia --project="$PKG" --machine-file=$(generate_hostfile) "$PKG/bin/membrane-analysis.jl" ...
# Existing worker processes are used instead of starting new ones.
