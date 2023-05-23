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
            tag_file = output_dir * "tag_1.h5"
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
        
        # calculating box dimensions and area

        box_dimensions(;
            traj_file=traj_file,
            area_file=output_dir * "A.dat",
            box_dims_file=output_dir * "box_dims.dat")
        
        @test readdlm(output_dir * "A.dat")[1] ≈ readdlm(traj_dir * "A.dat")[1]
        @test readdlm(output_dir * "box_dims.dat")[1:3] ≈ readdlm(traj_dir * "box_dims.dat")[1:3]

        # calculating |hq|^2 vs |q| data from fluctuation spectrum
        
        Lx, Ly = readdlm(output_dir * "box_dims.dat")[1:2]

        fluctuation_hq2_data(;
            box_dims=(Lx, Ly),
            fs_files=[fs_file],
            output_dir=output_dir,
            q_max=q_max)
        
        @test readdlm(output_dir * "hq2.dat")[:, 2] ≈ readdlm(traj_dir * "hq2.dat")[:, 2]
        
        # calculating blocking error
        
        @test blocking_error(sin.(1:1000)) ≈ 0.0012181125344636647
        
        # calculating index pairs
        
        @test get_index_pairs(2)[1] ≈ [0.0, 1.0, 1.4142135623730951, 2.0]
        @test get_index_pairs(2)[2] == [[(0, 0)], [(0, 1), (1, 0)], [(1, 1)], [(0, 2), (2, 0)]]
        
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

        # calculating reference atoms

        ref_atoms = find_ref_atoms(;
            input_dir=traj_dir,
            lipids=lipids,
            z=5,
        )

        @test ref_atoms == Dict(DOPE_m2 => "C3A", DOPC_m2 => "C3A")

        # calculating voronoi shells

        v_shells = voronoi_shells(;
           points=hcat([[i, j] for i in 1:2:9, j in 1:2:9]...),
           box_dims=(10, 10),
           center=[4, 4],
           n_shells=2
        )

        @test v_shells == Dict(2 => [1, 2, 3, 4, 6, 9, 11, 14, 16, 17, 18, 19], 1 => [7, 8, 12, 13])
        
        # calculating curvature spectrum of reference atoms of the lipids

        lipids_curvature_spectrum(;
            pdb_file=pdb_file,
            traj_files=[traj_file],
            fs_files=[fs_file],
            output_dir=output_dir,
            lipids=lipids,
            q_max=0.21,
            q_min=0.01,
            nqs=10)

        @test readdlm(output_dir * "DOPE_cqs.dat") ≈ readdlm(traj_dir * "DOPE_cqs.dat")
        @test readdlm(output_dir * "DOPC_cqs.dat") ≈ readdlm(traj_dir * "DOPC_cqs.dat")

        # calculating bending modulus from TCB method

        TCB_analysis(;
            input_dir=output_dir,
            lipids=lipids,
            z_cutoff=20,
            output_dir=output_dir,
        )

        @test readdlm(output_dir * "kc.dat")[1] ≈ readdlm(traj_dir * "kc.dat")[1]

        # calculating bending modulus from undulation spectrum method

        hq2_analysis(;
            input_dir=output_dir,
            output_dir=output_dir,
        )

        @test readdlm(output_dir * "kc_hq2.dat")[1] ≈ readdlm(traj_dir * "kc_hq2.dat")[1]

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

        # calculating RDF

        lipids_radial_distribution(;
            pdb_file=pdb_file,
            traj_file=traj_file,
            fs_file=fs_file,
            output_file=output_dir * "RDF.dat",
            lipids=lipids,
            lipid_A=DOPE_m2,
            lipid_B=DOPE_m2,
        )

        @test readdlm(output_dir * "RDF.dat") ≈ readdlm(traj_dir * "RDF.dat")
       
        # calculating number of H-bonds

        count_interactions(;
            pdb_file=pdb_file,
            traj_file=traj_file,
            fs_file=fs_file,
            output_file=tag_file,
            lipids=lipids,
        )
        
        @test isfile(tag_file)

        rm(output_dir, recursive=true)
end
