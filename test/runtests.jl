using MembraneAnalysis
using Test

@testset "MembraneAnalysis.jl" begin
        
        # setting system parameters

        begin
            const traj_dir = "./"
            const traj_file = traj_dir * "t_1.xtc"
            const pdb_file = traj_dir * "s_1.pdb"
            const lipids = [DOPE_m2, DOPC_m2]
            const L_grid = 15
            const q_max = 0.08
            const output_dir = "./temp/"
            if !isdir(output_dir) mkdir(output_dir) end
            const fs_file = output_dir * "fs_1.h5"
            const ds_file = output_dir * "ds_1.h5"
        end

        # calculating fluctuation spectrum of trajectories

        fluctuation_spectrum(;
            pdb_file=pdb_file,
            traj_file=traj_file,
            output_file=fs_file,
            lipids=lipids,
            L_grid=L_grid
        )

        # calculating area expansion modulus

        area_expansion_modulus(;
            traj_files=[traj_file],
            output_file=output_dir * "KA.dat"
        )

        # calculating mean height of lipids heavy atoms from the midplane

        lipids_atoms_height(;
            pdb_file=pdb_file,
            traj_file=traj_file,
            output_dir=output_dir,
            lipids=lipids
        )

        # calculating mean sampled curvature of heavy atoms of the lipids

        lipids_sampled_curvature(;
            pdb_file=pdb_file,
            traj_files=[traj_file],
            fs_files=[fs_file],
            output_dir=output_dir,
            lipids=lipids,
            q_max=q_max
        )

        # calculating curvature spectrum of reference atoms of the lipids

        lipids_curvature_spectrum(;
            pdb_file=pdb_file,
            traj_files=[traj_file],
            fs_files=[fs_file],
            output_dir=output_dir,
            lipids=lipids,
            q_max=q_max
        )

        # calculating thickness spectrum of reference atoms of the lipids

        lipids_thickness_spectrum(;
            pdb_file=pdb_file,
            traj_files=[traj_file],
            fs_files=[fs_file],
            output_dir=output_dir,
            lipids=lipids,
            q_max=q_max
        )

        # calculating density spectrum of trajectories

        lipids_density_spectrum(;
            pdb_file=pdb_file,
            traj_file=traj_file,
            output_file=ds_file,
            lipids=lipids,
            L_grid=L_grid
        )

        rm(output_dir)
        
end
