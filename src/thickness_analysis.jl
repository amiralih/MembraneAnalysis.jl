"""
    lipids_thickness_spectrum(;
        pdb_file,
        traj_files,
        fs_files,
        output_dir,
        lipids,
        ref_atoms=Dict(lipid => lipid.ref_atom for lipid in lipids),
        q_max)

Calculates thickness spectrum of the lipids using their reference atom position. Assumes square (Lx = Ly) bilayer. Results for lipid "XXXX" will be stored in `XXXX_tqs.dat` in the output directory.

### Keyword arguments

* `pdb_file`: PDB structure file;
* `traj_files`: a list of trajectory files;
* `fs_files`: a list of corresponding HDF5 fluctuation spectrum files of the trajectory files;
* `output_dir`: output directory;
* `lipids`: a list of lipids of type `Lipid` as defined in `lipids.jl`;
* `ref_atoms`: a dictionary of reference atoms for each lipid;
* `q_max`: maximum q mode magnitude value to be used.

"""
function lipids_thickness_spectrum(;
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
 
    # creating dictionary of vectors to record t_q values
    
    tqs = Dict()

    for lipid in lipids
        tqs[lipid.name] = Dict()
        for q_i in 2:length(qs)
            tqs[lipid.name][q_i] = Float64[]
        end
    end

    # reading trajectories    

    for (traj_file, fs_file) in zip(traj_files, fs_files)
        traj = Chemfiles.Trajectory(traj_file)
        n_frames = Int(Chemfiles.size(traj))
        
        tq = h5read(fs_file, "tq")
        l_id = h5read(fs_file, "l_id")

        println("Calculating thickness spectrum for trajectory file $(traj_file)")
        @showprogress 0.01 "Analyzing trajectory..." for frame_index in 1:n_frames
            
            # read a frame
            
            leaflet_id = l_id[frame_index, :]

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
                    total_t = 0.0
                   
                    # two degrees of freedom per mode

                    N_modes = 2 * length(index_pairs[q_i])
                    N_lipids = length(ref_inds) 

                    for index in ref_inds_1
                        X = coords[1, index]
                        Y = coords[2, index]
                        for (I_ind, J_ind) in index_pairs[q_i]
                            x_ind = I_ind + 1
                            y_ind = J_ind + 1
                            qx = I_ind * 2π / Lx
                            qy = J_ind * 2π / Ly
                            total_t += real((1 / (Lx * Ly)) * tq[frame_index, x_ind, y_ind] *
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
                            total_t += real((1 / (Lx * Ly)) * tq[frame_index, x_ind, y_ind] *
                                            exp(im * (qx * X + qy * Y))) 
                        end
                    end

                    push!(tqs[lipid.name][q_i], total_t / (N_modes * N_lipids))
                end
            end
        end
        Chemfiles.close(traj); GC.gc()
    end

    for lipid in lipids
        
        output_file = output_dir * lipid.name * "_tqs.dat"
        output_tqs = Float64[]
        output_tqes = Float64[]
        
        for q_i in 2:length(qs)
            push!(output_tqs, mean(tqs[lipid.name][q_i]))
            push!(output_tqes, blocking_error(tqs[lipid.name][q_i]))
        end
        
        writedlm(output_file, [qs[2:end] output_tqs output_tqes])
    end
    
    return nothing
end

"""
    peptide_thickness_spectrum(;
        pdb_file,
        traj_files,
        fs_files,
        output_file,
        lipids,
        ref_residue),
        q_max)

Calculates thickness spectrum of the peptide using the CA atom of its reference residue. Assumes square (Lx = Ly) bilayer.

### Keyword arguments

* `pdb_file`: PDB structure file;
* `traj_files`: a list of trajectory files;
* `fs_files`: a list of corresponding HDF5 fluctuation spectrum files of the trajectory files;
* `output_file`: output file;
* `lipids`: a list of lipids of type `Lipid` as defined in `lipids.jl`;
* `ref_residue`: residue number of the reference residue of the peptide;
* `q_max`: maximum q mode magnitude value to be used.

"""
function peptide_thickness_spectrum(;
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
 
    # creating dictionary of vectors to record t_q values
    
    tqs = Dict()

    tqs["PEP"] = Dict()

    for q_i in 2:length(qs)
        tqs["PEP"][q_i] = Float64[]
    end

    # reading trajectories    

    for (traj_file, fs_file) in zip(traj_files, fs_files)
        traj = Chemfiles.Trajectory(traj_file)
        n_frames = Int(Chemfiles.size(traj))
        
        tq = h5read(fs_file, "tq")
        l_id = h5read(fs_file, "l_id")

        println("Calculating thickness spectrum for trajectory file $(traj_file)")
        @showprogress 0.01 "Analyzing trajectory..." for frame_index in 1:n_frames
            
            # read a frame
            
            leaflet_id = l_id[frame_index, :]

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
                total_t = 0.0
            
                # two degrees of freedom per mode

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
                        total_t += real((1 / (Lx * Ly)) * tq[frame_index, x_ind, y_ind] *
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
                        total_t += real((1 / (Lx * Ly)) * tq[frame_index, x_ind, y_ind] *
                                        exp(im * (qx * X + qy * Y)))
                    end
                end
                push!(tqs["PEP"][q_i], total_t / (N_modes * N_peptides))
            end
        end
        Chemfiles.close(traj); GC.gc()
    end
    
    output_tqs = Float64[]
    output_tqes = Float64[]
    
    for q_i in 2:length(qs)
        push!(output_tqs, mean(tqs["PEP"][q_i]))
        push!(output_tqes, blocking_error(tqs["PEP"][q_i]))
    end
    
    writedlm(output_file, [qs[2:end] output_tqs output_tqes])
    
    return nothing
end

