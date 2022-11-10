"""
Calculates curvature per mode of a point (x, y) from fluctuation spectrum.
"""
function curvature(;
    point,
    hq,
    box_dims,
    q_max
)
    Lx, Ly = box_dims

    c = 0im
    
    # number of degrees of freedom (sines and cosines) = one per q value since q_(m=1, n=0) = [q_(m=-1, n=0)]^*
    
    n_dof = 0

    n_grid_x, n_grid_y = size(hq)
    n_max_x = Int(floor(n_grid_x / 2)) - 1
    n_max_y = Int(floor(n_grid_y / 2)) - 1
    
    for i in -n_max_x:n_max_x, j in -n_max_y:n_max_y
        qx, qy = (2π ./ (Lx, Ly)) .* (i, j)
        
        x_ind, y_ind = mod1.((i, j) .+ 1, (n_grid_x, n_grid_y))
        
        if (qx^2 + qy^2) < q_max^2
            c += (1 / (Lx * Ly)) * hq[x_ind, y_ind] * (qx^2 + qy^2) * exp(im * (qx * point[1] + qy * point[2]))
            n_dof += 1
        end
    end

    # number of modes (-1 to exclude (0, 0) from count) is half of the number of degrees of freedom
    
    n_modes = (n_dof - 1) / 2
    c_per_mode = real(c) / n_modes

    return  c_per_mode
end

"""
Calculates mean sampled curvature of C atoms of each lipid.
"""
function lipids_sampled_curvature(;
    pdb_file,
    traj_files,
    fs_files,
    output_dir,
    lipids,
    q_max
)

    # finding indices of tail atoms

    atoms = PDBTools.readPDB(pdb_file)
    leaflet_id = atom_leaflet(pdb_file=pdb_file, lipids=lipids)

    tail_atoms_inds = [a.index for a in atoms if (a.resname, a.name) in [(lipid.name, lipid.tail_atom) for lipid in lipids]]

    lipid_Cs = Dict()

    for lipid in lipids
        lipid_Cs[lipid.name] = unique([a.name for a in atoms if a.resname == lipid.name && a.name[1] in ['C', 'D', 'G']])
    end
    
    # creating dictionary of vectors to record <c> values
    
    cs = Dict()
    cs_frame = Dict()

    for lipid in lipids
        cs[lipid.name] = Dict()
        for C in lipid_Cs[lipid.name]
            cs[lipid.name][C] = Float64[]
        end
    end

    for lipid in lipids
        cs_frame[lipid.name] = Dict()
        for C in lipid_Cs[lipid.name]
            cs_frame[lipid.name][C] = Float64[]
        end
    end

    # reading trajectories    

    for (traj_file, fs_file) in zip(traj_files, fs_files)
        traj = Chemfiles.Trajectory(traj_file)
        n_frames = Int(Chemfiles.size(traj))
        
        hq = h5read(fs_file, "hq")

        println("Calculating mean curvatures for trajectory file $(traj_file)")
        @showprogress 0.01 "Analyzing trajectory..." for frame_index in 1:n_frames
            
            # read a frame

            frame = Chemfiles.read_step(traj, frame_index - 1)
            box_dims = Chemfiles.lengths(Chemfiles.UnitCell(frame))
            Lx, Ly, Lz = box_dims
            coords = Chemfiles.positions(frame)

            # making sure coordinates are within boundries
            
            coords = mod.(coords, box_dims)

            # fixing z values if tail atoms are not near the middle Z values

            while std(coords[3, tail_atoms_inds]) > 15
                coords[3, :] = mod.(coords[3, :] .+ 20, Lz)
            end

            midplane = mean(coords[3, tail_atoms_inds])
            
            Δz = (Lz / 2) - midplane
            if abs(Δz) > 0.05 * Lz
                coords[3, :] = mod.(coords[3, :] .+ Δz, Lz)
                midplane = mean(coords[3, tail_atoms_inds])
            end
            
            # calculating average curvature per mode of each atom

            for lipid in lipids
                for C in lipid_Cs[lipid.name]
                    cs_frame[lipid.name][C] = Float64[]
                end
            end

            for lipid in lipids
                for C in lipid_Cs[lipid.name]
                
                    # finding indices of carbon atoms of this lipid in each leaflet

                    C_inds = [a.index for a in atoms if a.resname == lipid.name && a.name == C] 
                    C_inds_1 = [i for i in C_inds if leaflet_id[i] == 1]
                    C_inds_2 = [i for i in C_inds if leaflet_id[i] == -1]

                    for index in C_inds_1
                        push!(cs_frame[lipid.name][C], curvature(point=coords[1:2, index],
                                                                 hq=hq[frame_index,:,:],
                                                                 box_dims=(Lx, Ly), q_max=q_max))
                    end
                    for index in C_inds_2
                        push!(cs_frame[lipid.name][C], -curvature(point=coords[1:2, index],
                                                                  hq=hq[frame_index,:,:],
                                                                  box_dims=(Lx, Ly), q_max=q_max))
                    end
                end
            end

            for lipid in lipids
                for C in lipid_Cs[lipid.name]
                    push!(cs[lipid.name][C], mean(cs_frame[lipid.name][C]))
                end
            end
        end

        Chemfiles.close(traj); GC.gc()
    end

    for lipid in lipids
        output_file = output_dir * lipid.name * "_cs.dat"
        output_cs = Float64[]
        output_ces = Float64[]
        for C in lipid_Cs[lipid.name]
            push!(output_cs, mean(cs[lipid.name][C]))
            push!(output_ces, blocking_error(cs[lipid.name][C]))
        end
        
        writedlm(output_file, [lipid_Cs[lipid.name] output_cs output_ces])
    end
    return nothing
end

"""
Calculate average curvature of CA atoms of peptide residues.
"""
function peptide_sampled_curvature(;
    pdb_file,
    traj_files,
    fs_files,
    output_file,
    lipids,
    n_residues,
    q_max
)
    # finding indices of tail atoms

    atoms = PDBTools.readPDB(pdb_file)
    leaflet_id = atom_leaflet(pdb_file=pdb_file, lipids=lipids)

    tail_atoms_inds = [a.index for a in atoms if (a.resname, a.name) in [(lipid.name, lipid.tail_atom) for lipid in lipids]]

    # creating dictionary of vectors to record <c> values
    
    cs = Dict()
    cs_frame = Dict()

    cs["PEP"] = Dict()

    for i in 1:n_residues
        cs["PEP"][i] = Float64[]
    end

    # reading trajectories    

    for (traj_file, fs_file) in zip(traj_files, fs_files)
        traj = Chemfiles.Trajectory(traj_file)
        n_frames = Int(Chemfiles.size(traj))
        
        hq = h5read(fs_file, "hq")

        println("Calculating mean curvature for trajectory file $(traj_file)")
        @showprogress 0.01 "Analyzing trajectory..." for frame_index in 1:n_frames
            
            # read a frame

            frame = Chemfiles.read_step(traj, frame_index - 1)
            box_dims = Chemfiles.lengths(Chemfiles.UnitCell(frame))
            Lx, Ly, Lz = box_dims
            coords = Chemfiles.positions(frame)

            # making sure coordinates are within boundries
            
            coords = mod.(coords, box_dims)

            # fixing z values if tail atoms are not near the middle Z values

            while std(coords[3, tail_atoms_inds]) > 15
                coords[3, :] = mod.(coords[3, :] .+ 20, Lz)
            end

            midplane = mean(coords[3, tail_atoms_inds])
            
            Δz = (Lz / 2) - midplane
            if abs(Δz) > 0.05 * Lz
                coords[3, :] = mod.(coords[3, :] .+ Δz, Lz)
                midplane = mean(coords[3, tail_atoms_inds])
            end
            
            # calculating average curvature per mode of each atom

            for i in 1:n_residues
            
                # finding indices of CA atoms of this residue in each leaflet

                inds = [a.index for a in atoms if a.segname[1:3] == "PRO" && a.resnum == i && a.name == "CA"]
                inds_1 = [j for j in inds if atoms[j].z > 0.0]
                inds_2 = [j for j in inds if atoms[j].z < 0.0]

                for index in inds_1
                    push!(cs["PEP"][i], curvature(point=coords[1:2, index],
                                                  hq=hq[frame_index,:,:],
                                                  box_dims=(Lx, Ly), q_max=q_max))
                end
                for index in inds_2
                    push!(cs["PEP"][i], -curvature(point=coords[1:2, index],
                                                   hq=hq[frame_index,:,:],
                                                   box_dims=(Lx, Ly), q_max=q_max))
                end
            end
        end

        Chemfiles.close(traj); GC.gc()
    end

    output_cs = Float64[]
    output_ces = Float64[]
    
    for i in 1:n_residues
        push!(output_cs, mean(cs["PEP"][i]))
        push!(output_ces, blocking_error(cs["PEP"][i]))
    end
    
    writedlm(output_file, [1:n_residues output_cs output_ces])
    
    return nothing
end

"""
Calculates curvature spectrum of the lipids using their reference atom position. Assumes square (Lx = Ly) bilayer.
"""
function lipids_curvature_spectrum(;
    pdb_file,
    traj_files,
    fs_files,
    output_dir,
    lipids,
    ref_atoms=Dict(lipid => lipid.ref_atom for lipid in lipids),
    q_max
)
    # finding indices of tail atoms

    atoms = PDBTools.readPDB(pdb_file)
    leaflet_id = atom_leaflet(pdb_file=pdb_file, lipids=lipids)

    tail_atoms_inds = [a.index for a in atoms if (a.resname, a.name) in [(lipid.name, lipid.tail_atom) for lipid in lipids]]

    # find |q| values (up to q_max) and corresponing indices
    
    (L, _, _) = box_dimensions(traj_file=traj_files[1])

    # including at least 9 q values (for N = 4) even if greater than q_max

    N = max(Int(floor(q_max / (2π / L))), 4)
    values, index_pairs = get_index_pairs(N)
    qs = values * (2π / L)
 
    # creating dictionary of vectors to record c_q values
    
    cqs = Dict()

    for lipid in lipids
        cqs[lipid.name] = Dict()
        for q_i in 2:length(qs)
            cqs[lipid.name][q_i] = Float64[]
        end
    end

    # reading trajectories    

    for (traj_file, fs_file) in zip(traj_files, fs_files)
        traj = Chemfiles.Trajectory(traj_file)
        n_frames = Int(Chemfiles.size(traj))
        
        hq = h5read(fs_file, "hq")

        println("Calculating curvature spectrum for trajectory file $(traj_file)")
        @showprogress 0.01 "Analyzing trajectory..." for frame_index in 1:n_frames
            
            # read a frame

            frame = Chemfiles.read_step(traj, frame_index - 1)
            box_dims = Chemfiles.lengths(Chemfiles.UnitCell(frame))
            Lx, Ly, Lz = box_dims
            coords = Chemfiles.positions(frame)

            # making sure coordinates are within boundries
            
            coords = mod.(coords, box_dims)

            # fixing z values if tail atoms are not near the middle Z values

            while std(coords[3, tail_atoms_inds]) > 15
                coords[3, :] = mod.(coords[3, :] .+ 20, Lz)
            end

            midplane = mean(coords[3, tail_atoms_inds])
            
            Δz = (Lz / 2) - midplane
            if abs(Δz) > 0.05 * Lz
                coords[3, :] = mod.(coords[3, :] .+ Δz, Lz)
                midplane = mean(coords[3, tail_atoms_inds])
            end
            
            # calculating average curvature per mode at |q| values for the reference atom of each lipid

            for lipid in lipids
            
                # finding indices of reference atoms of this lipid in each leaflet

                ref_inds = [a.index for a in atoms if a.resname == lipid.name && a.name == ref_atoms[lipid]] 
                ref_inds_1 = [i for i in ref_inds if leaflet_id[i] == 1]
                ref_inds_2 = [i for i in ref_inds if leaflet_id[i] == -1]
                
                for q_i in 2:length(qs)
                    total_c = 0.0
                    
                    # two degrees of freedom per mode

                    N_modes = length(index_pairs[q_i])
                    N_lipids = length(ref_inds) 

                    for index in ref_inds_1
                        X = coords[1, index]
                        Y = coords[2, index]
                        for (I_ind, J_ind) in index_pairs[q_i]
                            x_ind = I_ind + 1
                            y_ind = J_ind + 1
                            qx = I_ind * 2π / Lx
                            qy = J_ind * 2π / Ly
                            total_c += real((1 / (Lx * Ly)) * qs[q_i]^2 * hq[frame_index, x_ind, y_ind] *
                                            exp(im * (qx * X + qy * Y)))
                        end
                    end
                    for index in ref_inds_2
                        X = coords[1, index]
                        Y = coords[2, index]
                        for (I_ind, J_ind) in index_pairs[q_i]
                            x_ind = I_ind + 1
                            y_ind = J_ind + 1
                            qx = I_ind * 2π / Lx
                            qy = J_ind * 2π / Ly
                            total_c -= real((1 / (Lx * Ly)) * qs[q_i]^2 * hq[frame_index, x_ind, y_ind] *
                                            exp(im * (qx * X + qy * Y)))
                        end
                    end
            
                    push!(cqs[lipid.name][q_i], total_c / (N_modes * N_lipids))
                end
            end
        end
        Chemfiles.close(traj); GC.gc()
    end

    for lipid in lipids
        
        output_file = output_dir * lipid.name * "_cqs.dat"
        output_cqs = Float64[]
        output_cqes = Float64[]

        for q_i in 2:length(qs)
            push!(output_cqs, mean(cqs[lipid.name][q_i]))
            push!(output_cqes, blocking_error(cqs[lipid.name][q_i]))
        end
        
        writedlm(output_file, [qs[2:end] output_cqs output_cqes])
    end
    
    return nothing
end

"""
Calculates curvature spectrum of the peptide using the CA atom of its reference residue. Assumes square (Lx = Ly) bilayer.
"""
function peptide_curvature_spectrum(;
    pdb_file,
    traj_files,
    fs_files,
    output_file,
    lipids,
    ref_residue,
    q_max
)
    # finding indices of tail atoms

    atoms = PDBTools.readPDB(pdb_file)
    leaflet_id = atom_leaflet(pdb_file=pdb_file, lipids=lipids)

    tail_atoms_inds = [a.index for a in atoms if (a.resname, a.name) in [(lipid.name, lipid.tail_atom) for lipid in lipids]]

    # find |q| values (up to q_max) and corresponing indices
    
    (L, _, _) = box_dimensions(traj_file=traj_files[1])

    # including at least 9 q values (for N = 4) even if greater than q_max

    N = max(Int(floor(q_max / (2π / L))), 4)
    values, index_pairs = get_index_pairs(N)
    qs = values * (2π / L)
 
    # creating dictionary of vectors to record c_q values
    
    cqs = Dict()

    cqs["PEP"] = Dict()

    for q_i in 2:length(qs)
        cqs["PEP"][q_i] = Float64[]
    end

    # reading trajectories    

    for (traj_file, fs_file) in zip(traj_files, fs_files)
        traj = Chemfiles.Trajectory(traj_file)
        n_frames = Int(Chemfiles.size(traj))
        
        hq = h5read(fs_file, "hq")

        println("Calculating curvature spectrum for trajectory file $(traj_file)")
        @showprogress 0.01 "Analyzing trajectory..." for frame_index in 1:n_frames
            
            # read a frame

            frame = Chemfiles.read_step(traj, frame_index - 1)
            box_dims = Chemfiles.lengths(Chemfiles.UnitCell(frame))
            Lx, Ly, Lz = box_dims
            coords = Chemfiles.positions(frame)

            # making sure coordinates are within boundries
            
            coords = mod.(coords, box_dims)

            # fixing z values if tail atoms are not near the middle Z values

            while std(coords[3, tail_atoms_inds]) > 15
                coords[3, :] = mod.(coords[3, :] .+ 20, Lz)
            end

            midplane = mean(coords[3, tail_atoms_inds])
            
            Δz = (Lz / 2) - midplane
            if abs(Δz) > 0.05 * Lz
                coords[3, :] = mod.(coords[3, :] .+ Δz, Lz)
                midplane = mean(coords[3, tail_atoms_inds])
            end
            
            # finding indices of CA atoms of the reference residue in each leaflet

            ref_inds = [a.index for a in atoms if a.segname[1:3] == "PRO" && a.resnum == ref_residue && a.name == "CA"]
            ref_inds_1 = [j for j in ref_inds if atoms[j].z > 0.0]
            ref_inds_2 = [j for j in ref_inds if atoms[j].z < 0.0]

            for q_i in 2:length(qs)
                total_c = 0.0
                
                # tow degrees of freedom per mode

                N_modes = 2 * length(index_pairs[q_i])
                N_peptides = length(ref_inds) 

                for index in ref_inds_1
                    X = coords[1, index]
                    Y = coords[2, index]
                    for (I_ind, J_ind) in index_pairs[q_i]
                        x_ind = I_ind + 1
                        y_ind = J_ind + 1
                        qx = I_ind * 2π / Lx
                        qy = J_ind * 2π / Ly
                        total_c += real((1 / (Lx * Ly)) * qs[q_i]^2 * hq[frame_index, x_ind, y_ind] *
                                        exp(im * (qx * X + qy * Y)))
                    end
                end
                for index in ref_inds_2
                    X = coords[1, index]
                    Y = coords[2, index]
                    for (I_ind, J_ind) in index_pairs[q_i]
                        x_ind = I_ind + 1
                        y_ind = J_ind + 1
                        qx = I_ind * 2π / Lx
                        qy = J_ind * 2π / Ly
                        total_c -= real((1 / (Lx * Ly)) * qs[q_i]^2 * hq[frame_index, x_ind, y_ind] *
                                        exp(im * (qx * X + qy * Y)))
                    end
                end
                push!(cqs["PEP"][q_i], total_c / (N_modes * N_peptides))
            end
        end
        Chemfiles.close(traj); GC.gc()
    end
    
    output_cqs = Float64[]
    output_cqes = Float64[]
    
    for q_i in 2:length(qs)
        push!(output_cqs, mean(cqs["PEP"][q_i]))
        push!(output_cqes, blocking_error(cqs["PEP"][q_i]))
    end
    
    writedlm(output_file, [qs[2:end] output_cqs output_cqes])
    
    return nothing
end

