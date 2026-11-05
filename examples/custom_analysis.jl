# Using MembraneAnalysis.jl from Julia: running the pipeline and writing a custom analysis.
#
#   julia --project=/path/to/MembraneAnalysis.jl custom_analysis.jl

using MembraneAnalysis

pdb = "step5_charmm2gmx.pdb"
trajectories = ["step7_$(i).xtc" for i in 1:15]

# 1. The same as the command-line tool
summary = run_analysis("parameters.toml", pdb, trajectories; output_dir="results", nprocs=8)
println("K_A = ", summary["box"]["KA"], " ± ", summary["box"]["KA_error"], " kBT/Å²")

# 2. A custom per-frame quantity with the memory-safe frame loop: the curvature sampled by
#    the backbone bead of peptide residue 4 in every frame (needs the spectra, so surface
#    settings are passed).
lipids = [Lipid("POPC", "NC3", "C1A", "C4B", 12)]
proteins = [Protein("PEP"; segments=["PROA"])]
system = MembraneSystem(pdb, lipids, proteins)
surface = SurfaceSettings(grid_size=(14, 14), fill=:min_curvature,
                          mask_atoms=system.protein_backbone[1], mask_radius=5.0)

bb4 = MembraneAnalysis.protein_residue_group(system, 1, 4).atoms[1]
leaflet_sign = system.protein_sign[bb4]  # +1 if the peptide is in the upper leaflet
curvatures = Float64[]
foreach_frame(system, trajectories[1]; surface=surface) do ctx, i
    # curvature per mode from modes with 0 < |q| < 0.08 1/Å
    x, y = ctx.coords[1, bb4], ctx.coords[2, bb4]
    push!(curvatures, leaflet_sign * curvature_at(ctx.hq, ctx.Lx, ctx.Ly, x, y; q_max=0.08))
end
println("mean curvature at residue 4: ", sum(curvatures) / length(curvatures),
        " ± ", blocking_error(curvatures), " 1/Å")
