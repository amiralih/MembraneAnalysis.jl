using MembraneAnalysis
using Test
using DelimitedFiles

@testset "MembraneAnalysis.jl" begin
        
        # setting system parameters

        begin
            traj_dir = "./sample_files/"
            traj_file = traj_dir * "t_1.xtc"
            pdb_file = traj_dir * "s_1.pdb"
            lipids = [DOPE_m2, DOPC_m2]
            L_grid = 15
            q_max = 0.08
            output_dir = "./temp_out/"
            if !isdir(output_dir) mkdir(output_dir) end
            fs_file = output_dir * "fs_1.h5"
            ds_file = output_dir * "ds_1.h5"
        end

        # calculating fluctuation spectrum of trajectories

        fluctuation_spectrum(;
            pdb_file=pdb_file,
            traj_file=traj_file,
            output_file=fs_file,
            lipids=lipids,
            L_grid=L_grid
        )

        @test isfile(fs_file)

        # calculating area expansion modulus

        area_expansion_modulus(;
            traj_files=[traj_file],
            output_file=output_dir * "KA.dat"
        )

        @test readdlm(output_dir * "KA.dat") ≈ readdlm(traj_dir * "KA.dat")

        # calculating mean height of lipids heavy atoms from the midplane

        lipids_atoms_height(;
            pdb_file=pdb_file,
            traj_file=traj_file,
            fs_file=fs_file,
            output_dir=output_dir,
            lipids=lipids
        )

        @test readdlm(output_dir * "DOPE_zs.dat")[:, 2] ≈ readdlm(traj_dir * "DOPE_zs.dat")[:, 2]
        @test readdlm(output_dir * "DOPC_zs.dat")[:, 2] ≈ readdlm(traj_dir * "DOPC_zs.dat")[:, 2]
        
        # calculating mean sampled curvature of heavy atoms of the lipids

        lipids_sampled_curvature(;
            pdb_file=pdb_file,
            traj_files=[traj_file],
            fs_files=[fs_file],
            output_dir=output_dir,
            lipids=lipids,
            q_max=q_max
        )

        @test readdlm(output_dir * "DOPE_cs.dat")[:, 2] ≈ readdlm(traj_dir * "DOPE_cs.dat")[:, 2]
        @test readdlm(output_dir * "DOPC_cs.dat")[:, 2] ≈ readdlm(traj_dir * "DOPC_cs.dat")[:, 2]
#=        
        # calculating curvature spectrum of reference atoms of the lipids

        lipids_curvature_spectrum(;
            pdb_file=pdb_file,
            traj_files=[traj_file],
            fs_files=[fs_file],
            output_dir=output_dir,
            lipids=lipids,
            q_max=q_max
        )

        @test readdlm(output_dir * "DOPE_cqs.dat") ≈ readdlm(traj_dir * "DOPE_cqs.dat")
        @test readdlm(output_dir * "DOPC_cqs.dat") ≈ readdlm(traj_dir * "DOPC_cqs.dat")
=#        
        # calculating thickness spectrum of reference atoms of the lipids

        lipids_thickness_spectrum(;
            pdb_file=pdb_file,
            traj_files=[traj_file],
            fs_files=[fs_file],
            output_dir=output_dir,
            lipids=lipids,
            q_max=q_max
        )

        @test readdlm(output_dir * "DOPE_tqs.dat") ≈ readdlm(traj_dir * "DOPE_tqs.dat")
        @test readdlm(output_dir * "DOPC_tqs.dat") ≈ readdlm(traj_dir * "DOPC_tqs.dat")
        
        # calculating density spectrum of trajectories

        lipids_density_spectrum(;
            pdb_file=pdb_file,
            traj_file=traj_file,
            fs_file=fs_file,
            output_file=ds_file,
            lipids=lipids,
            L_grid=L_grid
        )

        @test isfile(ds_file)

        rm(output_dir, recursive=true)

end
