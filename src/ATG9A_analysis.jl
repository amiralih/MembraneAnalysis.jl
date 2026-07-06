function membrane_sampled_height(;
    pdb_file,
    traj_files,
    output_dir,
    lipids,
    ref_top,
    ref_left,
    ref_right,
    ref_center
)

    # finding indices of tail atoms

    atoms = PDBTools.readPDB(pdb_file)
    leaflet_id = atom_leaflet(pdb_file=pdb_file, lipids=lipids)

    tail_atoms_inds = [a.index for a in atoms if (a.resname, a.name) in [(lipid.name, lipid.tail_atom) for lipid in lipids]]
    ref_atoms_inds = [a.index for a in atoms if (a.resname, a.name) in [(lipid.name, "PO4") for lipid in lipids]]
    
    ref_ind_top = [a.index for a in atoms if a.segname[1:3] == "PRO" && a.residue == ref_top && a.name == "BB"][1]
    ref_ind_left = [a.index for a in atoms if a.segname[1:3] == "PRO" && a.residue == ref_left && a.name == "BB"][1]
    ref_ind_right = [a.index for a in atoms if a.segname[1:3] == "PRO" && a.residue == ref_right && a.name == "BB"][1]
    ref_ind_center = [a.index for a in atoms if a.segname[1:3] == "PRO" && a.residue == ref_center && a.name == "BB"][1]
    
    # store lipid height and number data
    
    n_1 = zeros(Int, 20, 20)
    z_1 = zeros(Float32, 20, 20)
    n_2 = zeros(Int, 20, 20)
    z_2 = zeros(Float32, 20, 20)
     
    # reading trajectories    

    for traj_file in traj_files
        traj = Chemfiles.Trajectory(traj_file)
        n_frames = Int(Chemfiles.size(traj))
        
        println("Processing trajectory file $(traj_file)")
        @showprogress 0.01 "Analyzing trajectory..." for frame_index in 1:n_frames
            
            # read a frame
            
            frame = Chemfiles.read_step(traj, frame_index - 1)
            box_dims = Chemfiles.lengths(Chemfiles.UnitCell(frame))
            Lx, Ly, Lz = box_dims
            coords = Chemfiles.positions(frame)

            # making sure coordinates are within boundries
            
            coords = mod.(coords, box_dims)

            # fixing z values if tail atoms are not near the middle Z values

            while std(coords[3, tail_atoms_inds]) > 20
                coords[3, :] = mod.(coords[3, :] .+ 20, Lz)
            end

            midplane = mean(coords[3, tail_atoms_inds])
            
            Δz = (Lz / 2) - midplane
            if abs(Δz) > 0.05 * Lz
                coords[3, :] = mod.(coords[3, :] .+ Δz, Lz)
                midplane = mean(coords[3, tail_atoms_inds])
            end

            m_x = coords[1, ref_ind_center]
            Δx = (Lx / 2) - m_x
            if abs(Δx) > 0.05 * Lx
                coords[1, :] = mod.(coords[1, :] .+ Δx, Lx)
            end

            m_y = coords[2, ref_ind_center]
            Δy = (Ly / 2) - m_y
            if abs(Δy) > 0.05 * Ly
                coords[2, :] = mod.(coords[2, :] .+ Δy, Ly)
            end
            
            m_x = (coords[1, ref_ind_top] + coords[1, ref_ind_left] + coords[1, ref_ind_right]) / 3
            m_y = (coords[2, ref_ind_top] + coords[2, ref_ind_left] + coords[2, ref_ind_right]) / 3
            
            leaflet_id = zeros(Int, length(atoms))
            
            for lipid in lipids
                indices = [a.index for a in atoms if a.resname == lipid.name && a.name == lipid.head_atom]
                for ind in indices
                    if coords[3, ind] > coords[3, ind + lipid.n_atoms - 1]
                        leaflet_id[ind:(ind + lipid.n_atoms - 1)] .= 1
                    else
                        leaflet_id[ind:(ind + lipid.n_atoms - 1)] .= -1
                    end
                end
            end
            
            # calculating polar bin of each lipid

            for ind in ref_atoms_inds
                a_x = coords[1, ind]
                a_y = coords[2, ind]
                a_z = coords[3, ind]
            
                t_x = coords[1, ref_ind_top]
                t_y = coords[2, ref_ind_top]
                
                a_rad = peuclidean((a_x, a_y), (m_x, m_y), (Lx, Ly))
                if 30.0 < a_rad < 130.0
                    b_rad = Int(floor((a_rad - 30.0) / 5.0) + 1)
                    
                    a_ang = atand(a_y - m_y, a_x - m_x)
                    t_ang = atand(t_y - m_y, t_x - m_x)
                    b_ang = a_ang - t_ang
                    if b_ang < 0.0
                        b_ang += 360.0
                    end
                    b_ang = Int(floor((b_ang / 360.0) * 20.0) + 1)
                    
                    if leaflet_id[ind] > 0
                        n_1[b_rad, b_ang] += 1
                        z_1[b_rad, b_ang] += a_z - midplane
                    end
                    if leaflet_id[ind] < 0
                        n_2[b_rad, b_ang] += 1
                        z_2[b_rad, b_ang] += a_z - midplane
                    end
                end
            end
        end

        Chemfiles.close(traj); GC.gc()
    end
     
    output_file = output_dir * "h_data.h5"
    if ispath(output_file) rm(output_file) end 
    h5write(output_file, "n_1", n_1)
    h5write(output_file, "z_1", z_1)
    h5write(output_file, "n_2", n_2)
    h5write(output_file, "z_2", z_2)
    
    return nothing
end

