"""
    thickness(;
        point,
        tq,
        box_dims,
        q_max,
        q_min=0)

Calculates thickness per mode of a point from fluctuation spectrum.

### Keyword arguments

* `point`: ordered pair of X and Y values;
* `tq`: 2D matrix of thickness fluctuation spectrum;
* `box_dims`: ordered pair of simulation box Lx and Ly values;
* `q_max`: maximum q mode magnitude value to be used.
* `q_min`: minimum q mode magnitude value to be used (zero by default).

"""
function thickness(;
    point,
    tq,
    box_dims,
    q_max,
    q_min=0
)
    Lx, Ly = box_dims

    t = 0im
    
    # number of degrees of freedom (sines and cosines) = one per q value since q_(m=1, n=0) = [q_(m=-1, n=0)]^*
    
    n_dof = 0

    n_grid_x, n_grid_y = size(tq)
    n_max_x = Int(floor(n_grid_x / 2)) - 1
    n_max_y = Int(floor(n_grid_y / 2)) - 1
    
    for i in -n_max_x:n_max_x, j in -n_max_y:n_max_y
        qx, qy = (2π ./ (Lx, Ly)) .* (i, j)
        
        x_ind, y_ind = mod1.((i, j) .+ 1, (n_grid_x, n_grid_y))
        
        if q_min^2 < (qx^2 + qy^2) < q_max^2
            t += (1 / (Lx * Ly)) * tq[x_ind, y_ind] * exp(im * (qx * point[1] + qy * point[2]))
            n_dof += 1
        end
    end

    # number of modes is half of the number of degrees of freedom
    
    n_modes = n_dof / 2
    t_per_mode = real(t) / n_modes

    return  t_per_mode
end

"""
    lipids_thickness_spectrum(;
        pdb_file,
        traj_files,
        fs_files,
        output_dir,
        lipids,
        ref_atoms=Dict(lipid => lipid.ref_atom for lipid in lipids),
        q_max,
        q_min,
        nqs)

Calculates thickness spectrum of the lipids using their reference atom position. Assumes square (Lx = Ly) bilayer. Results for lipid "XXXX" will be stored in `XXXX_tqs.dat` in the output directory.

### Keyword arguments

* `pdb_file`: PDB structure file;
* `traj_files`: a list of trajectory files;
* `fs_files`: a list of corresponding HDF5 fluctuation spectrum files of the trajectory files;
* `output_dir`: output directory;
* `lipids`: a list of lipids of type `Lipid` as defined in `lipids.jl`;
* `ref_atoms`: a dictionary of reference atoms for each lipid;
* `q_max`: maximum q mode magnitude value to be used;
* `q_min`: minimum q mode magnitude value to be used;
* `nqs`: number of segments to divide the q range from q_min to q_max;
* `tag_files`: a list of corresponding lipid interaction tag files of the trajectory files. Uses all lipids by default;
* `tag_id`: the tag id of the selected lipids.

"""
function lipids_thickness_spectrum(;
    pdb_file,
    traj_files,
    fs_files,
    output_dir,
    lipids,
    ref_atoms=Dict(lipid => lipid.ref_atom for lipid in lipids),
    q_max=0.21,
    q_min=0.01,
    nqs=10,
    tag_files=false,
    tag_id=0
)
    # finding indices of tail atoms

    atoms = PDBTools.readPDB(pdb_file)

    tail_atoms_inds = [a.index for a in atoms if (a.resname, a.name) in [(lipid.name, lipid.tail_atom) for lipid in lipids]]

    # find mid-segment |q| values (q_min to q_max)

    dq = (q_max - q_min) / nqs
    qs = collect((q_min + (dq / 2)):dq:q_max)
    
    # creating dictionary of vectors to record t_q values
    
    tqs = Dict()

    for lipid in lipids
        tqs[lipid.name] = Dict()
        for q_i in 1:nqs
            tqs[lipid.name][q_i] = Float64[]
        end
    end
    
    # defining a flag for if there are tag files

    tag_flag = true
    if tag_files == false
        tag_flag = false
        tag_files = ["---" for fs_file in fs_files]
    end

    # create index to resnum dictionary

    index_to_resnum = Dict()

    for a in atoms
        if a.resname in [lipid.name for lipid in lipids]
            index_to_resnum[a.index] = a.resnum
        end
    end

    # reading trajectories    

    for (traj_file, fs_file, tag_file) in zip(traj_files, fs_files, tag_files)
        traj = Chemfiles.Trajectory(traj_file)
        n_frames = Int(Chemfiles.size(traj))
        
        tq = h5read(fs_file, "tq")
        l_id = h5read(fs_file, "l_id")
        if tag_flag tag = h5read(tag_file, "tag") end

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

            while std(coords[3, tail_atoms_inds]) > 20
                coords[3, :] = mod.(coords[3, :] .+ 20, Lz)
            end

            midplane = mean(coords[3, tail_atoms_inds])
            
            Δz = (Lz / 2) - midplane
            if abs(Δz) > 0.05 * Lz
                coords[3, :] = mod.(coords[3, :] .+ Δz, Lz)
                midplane = mean(coords[3, tail_atoms_inds])
            end
            
            # calculating average thickness per mode at |q| values for the reference atom of each lipid

            for lipid in lipids
            
                # finding indices of reference atoms of this lipid in each leaflet

                ref_inds = [a.index for a in atoms if a.resname == lipid.name && a.name == ref_atoms[lipid]] 
                
                # selcting the tagged lipids

                if tag_flag
                    ref_inds = [i for i in ref_inds if tag[frame_index, index_to_resnum[i]] == tag_id]
                end

                ref_inds_1 = [i for i in ref_inds if leaflet_id[i] == 1]
                ref_inds_2 = [i for i in ref_inds if leaflet_id[i] == -1]
                
                for q_i in 1:nqs

                    q2 = q_min + q_i * dq 
                    q1 = q_min + (q_i - 1) * dq

                    for index in ref_inds_1
                        push!(tqs[lipid.name][q_i], thickness(point=coords[1:2, index],
                                                                 tq=tq[frame_index,:,:],
                                                                 box_dims=(Lx, Ly),
                                                                 q_max=q2, q_min=q1))
                    end
                    
                    for index in ref_inds_2
                        push!(tqs[lipid.name][q_i], thickness(point=coords[1:2, index],
                                                                 tq=tq[frame_index,:,:],
                                                                 box_dims=(Lx, Ly),
                                                                 q_max=q2, q_min=q1))
                    end
                end
            end
        end
        Chemfiles.close(traj); GC.gc()
    end

    for lipid in lipids
        
        output_file = output_dir * lipid.name * "_tqs.dat"
        output_tqs = Float64[]
        output_tqes = Float64[]

        for q_i in 1:nqs
            push!(output_tqs, mean(tqs[lipid.name][q_i]))
            push!(output_tqes, blocking_error(tqs[lipid.name][q_i]))
        end
        
        writedlm(output_file, [qs output_tqs output_tqes])
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
        q_max,
        q_min,
        nqs
)

Calculates thickness spectrum of the peptide using the CA atom of its reference residue. Assumes square (Lx = Ly) bilayer.

### Keyword arguments

* `pdb_file`: PDB structure file;
* `traj_files`: a list of trajectory files;
* `fs_files`: a list of corresponding HDF5 fluctuation spectrum files of the trajectory files;
* `output_file`: output file;
* `lipids`: a list of lipids of type `Lipid` as defined in `lipids.jl`;
* `ref_residue`: residue number of the reference residue of the peptide;
* `q_max`: maximum q mode magnitude value to be used;
* `q_min`: minimum q mode magnitude value to be used;
* `nqs`: number of segments to divide the q range from q_min to q_max.

"""
function peptide_thickness_spectrum(;
    pdb_file,
    traj_files,
    fs_files,
    output_file,
    lipids,
    ref_residue,
    q_max=0.21,
    q_min=0.01,
    nqs=10
)
    # finding indices of tail atoms

    atoms = PDBTools.readPDB(pdb_file)

    tail_atoms_inds = [a.index for a in atoms if (a.resname, a.name) in [(lipid.name, lipid.tail_atom) for lipid in lipids]]
    
    # find mid-segment |q| values (q_min to q_max)

    dq = (q_max - q_min) / nqs
    qs = collect((q_min + (dq / 2)):dq:q_max)

    # creating dictionary of vectors to record c_q values
    
    tqs = Dict()

    for q_i in 1:nqs
        tqs[q_i] = Float64[]
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

            while std(coords[3, tail_atoms_inds]) > 20
                coords[3, :] = mod.(coords[3, :] .+ 20, Lz)
            end

            midplane = mean(coords[3, tail_atoms_inds])
            
            Δz = (Lz / 2) - midplane
            if abs(Δz) > 0.05 * Lz
                coords[3, :] = mod.(coords[3, :] .+ Δz, Lz)
                midplane = mean(coords[3, tail_atoms_inds])
            end
            
            # finding indices of CA atoms of the reference residue in each leaflet

            ref_inds = [a.index for a in atoms if a.segname[1:3] == "PRO" && a.residue == ref_residue && (a.name == "BB" || a.name == "CA")]
            ref_inds_1 = [j for j in ref_inds if atoms[j].z > 0.0]
            ref_inds_2 = [j for j in ref_inds if atoms[j].z < 0.0]

            for q_i in 1:nqs
                q2 = q_min + q_i * dq 
                q1 = q_min + (q_i - 1) * dq

                for index in ref_inds_1
                    push!(tqs[q_i], thickness(point=coords[1:2, index],
                                              tq=tq[frame_index,:,:],
                                              box_dims=(Lx, Ly),
                                              q_max=q2, q_min=q1))
                end
                
                for index in ref_inds_2
                    push!(tqs[q_i], thickness(point=coords[1:2, index],
                                              tq=tq[frame_index,:,:],
                                              box_dims=(Lx, Ly),
                                              q_max=q2, q_min=q1))
                end
            end
        end
        Chemfiles.close(traj); GC.gc()
    end
    
    output_tqs = Float64[]
    output_tqes = Float64[]
    
    for q_i in 1:nqs
        push!(output_tqs, mean(tqs[q_i]))
        push!(output_tqes, blocking_error(tqs[q_i]))
    end
    
    writedlm(output_file, [qs output_tqs output_tqes])
    
    return nothing
end

