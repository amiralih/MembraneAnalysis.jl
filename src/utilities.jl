"""
Determines the leaflet of each atom of the lipid in a bilayer based on vertical position of the head atoms in the PDB structure file. Outputs an array of ±1s.
"""
function atom_leaflet(;
        pdb_file,
        lipids,
)
    atoms = PDBTools.readPDB(pdb_file)
    
    leaflet_id = zeros(Int, length(atoms))
    
    for lipid in lipids
        indices = [a.index for a in atoms if a.resname == lipid.name && a.name == lipid.head_atom]
        for ind in indices
            if atoms[ind].z > atoms[ind + lipid.n_atoms - 1].z
                leaflet_id[ind:(ind + lipid.n_atoms - 1)] .= 1
            else
                leaflet_id[ind:(ind + lipid.n_atoms - 1)] .= -1
            end
        end
    end

    return leaflet_id
end

"""
Determines the leaflet of each atom of the lipid in a bilayer based on vertical position of the head atoms relative to the leaflet mean height. Outputs an array of ±1s.
"""
function atom_leaflet_dynamic(;
        coords,
        zs_m,
        Ls,
        pdb_file,
        lipids,
)
    atoms = PDBTools.readPDB(pdb_file)
    
    leaflet_id = zeros(Int, length(atoms))
    
    n_grid_x, n_grid_y = size(zs_m) 
    Lx, Ly = Ls

    for lipid in lipids
        indices = [a.index for a in atoms if a.resname == lipid.name && a.name == lipid.head_atom]
        for ind in indices
            x_index = Int(floor(coords[1, ind] * n_grid_x / Lx)) + 1
            y_index = Int(floor(coords[2, ind] * n_grid_y / Ly)) + 1
            if atoms[ind].z > zs_m[x_index, y_index]
                leaflet_id[ind:(ind + lipid.n_atoms - 1)] .= 1
            else
                leaflet_id[ind:(ind + lipid.n_atoms - 1)] .= -1
            end
        end
    end

    return leaflet_id
end

"""
Finds the mean dimensions of the simulation box.
"""
function box_dimensions(; traj_file)
    traj = Chemfiles.Trajectory(traj_file)
    n_frames = Int(Chemfiles.size(traj))
    
    Lxs = Float64[]
    Lys = Float64[]
    Lzs = Float64[]

    println("Finding simulation box dimensions for trajectory file $(traj_file)")
    @showprogress 0.01 "Analyzing trajectory..." for frame_index in 1:n_frames
        frame = Chemfiles.read_step(traj, frame_index - 1)
        box_dims = Chemfiles.lengths(Chemfiles.UnitCell(frame))
        push!(Lxs, box_dims[1])
        push!(Lys, box_dims[2])
        push!(Lzs, box_dims[3])
    end
    
    Chemfiles.close(traj); GC.gc()
    
    Lx = mean(Lxs)
    Ly = mean(Lys)
    Lz = mean(Lzs)

    return (Lx, Ly, Lz)
end

"""
Finds the standard error of the mean for an array divided in n blocks.
"""
function blocking_error(a, n=10)
    n_data = length(a)
    n_in_block = floor(Int, n_data / n)
    
    if n_in_block < 1
        error("too few ($(n_data)) data points for $(n) blocks")
    else
        a_blocked = zeros(eltype(a), n)
        
        for i in 1:n
            a_blocked[i] = mean(a[(1 + (i - 1) * n_in_block):(i * n_in_block)])
        end

        return std(a_blocked) / √n
    end
end

"""
Calculates the average distance of each carbon atom of each lipid from the midplane.
"""
function lipids_atoms_height(;
    pdb_file,
    traj_file,
    fs_file,
    output_dir,
    lipids
)

    # finding indices of tail atoms and carbons

    atoms = PDBTools.readPDB(pdb_file)
    leaflet_id = atom_leaflet(pdb_file=pdb_file, lipids=lipids)
    l_id = h5read(fs_file, "l_id")

    tail_atoms_inds = [a.index for a in atoms if (a.resname, a.name) in [(lipid.name, lipid.tail_atom) for lipid in lipids]]

    lipid_Cs = Dict()

    for lipid in lipids
        lipid_Cs[lipid.name] = unique([a.name for a in atoms if a.resname == lipid.name && a.name[1] in ['C', 'D', 'G']])
    end
    
    # creating dictionary of vectors to record z values
    
    zs = Dict()

    for lipid in lipids
        zs[lipid.name] = Dict()
        for C in lipid_Cs[lipid.name]
            zs[lipid.name][C] = Float64[]
        end
    end

    # opening trajectory

    traj = Chemfiles.Trajectory(traj_file)
    n_frames = Int(Chemfiles.size(traj))

    println("Finding z of lipid atoms for trajectory file $(traj_file)")
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
        
        # collecting z values

        for lipid in lipids
            for C in lipid_Cs[lipid.name]
            
                # finding indices of carbon atoms of this lipid in each leaflet

                C_inds = [a.index for a in atoms if a.resname == lipid.name && a.name == C] 
                C_inds_1 = [i for i in C_inds if leaflet_id[i] == 1]
                C_inds_2 = [i for i in C_inds if leaflet_id[i] == -1]

                append!(zs[lipid.name][C], coords[3, C_inds_1] .- midplane)
                append!(zs[lipid.name][C], midplane .- coords[3, C_inds_2])
            end
        end
    end

    Chemfiles.close(traj); GC.gc()

    for lipid in lipids
        output_file = output_dir * lipid.name * "_zs.dat"
        output_zs = Float64[]
        output_zes = Float64[]
        for C in lipid_Cs[lipid.name]
            push!(output_zs, mean(zs[lipid.name][C]))
            push!(output_zes, blocking_error(zs[lipid.name][C]))
        end
        
        writedlm(output_file, [lipid_Cs[lipid.name] output_zs output_zes])
    end
    
    return nothing
end

"""
Calculate average distance of CA atoms of each residue of the peptide from midplane.
"""
function peptide_atoms_height(;
    pdb_file,
    traj_files,
    output_file,
    lipids,
    n_residues
)

    # finding indices tail atoms

    atoms = PDBTools.readPDB(pdb_file)

    tail_atoms_inds = [a.index for a in atoms if (a.resname, a.name) in [(lipid.name, lipid.tail_atom) for lipid in lipids]]
    
    # creating dictionary of vectors to record z values of peptide amino acid residues
    
    zs = Dict()

    for i in 1:n_residues
        zs[i] = Float64[]
    end
    
    for traj_file in traj_files
            
        # opening trajectory

        traj = Chemfiles.Trajectory(traj_file)
        n_frames = Int(Chemfiles.size(traj))

        println("Finding z of peptide atoms for trajectory file $(traj_file)")
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
            
            # collecting z values

            for i in 1:n_residues
            
                # finding indices of CA atoms of this residue in each leaflet

                inds = [a.index for a in atoms if a.segname[1:3] == "PRO" && a.resnum == i && a.name == "CA"]
                inds_1 = [j for j in inds if atoms[j].z > 0.0]
                inds_2 = [j for j in inds if atoms[j].z < 0.0]

                append!(zs[i], coords[3, inds_1] .- midplane)
                append!(zs[i], midplane .- coords[3, inds_2])
            end
        end

        Chemfiles.close(traj); GC.gc()
    end

    output_zs = Float64[]
    output_zes = Float64[]
    
    for i in 1:n_residues
        push!(output_zs, mean(zs[i]))
        push!(output_zes, blocking_error(zs[i]))
    end
    
    writedlm(output_file, [1:n_residues output_zs output_zes])
    
    return nothing
end

"""
Find all (n, m) index pairs and corresponding v = √(n^2 + m^2) values where v ≤ N.
"""
function get_index_pairs(N)
    all_values = []
    all_index_pairs = []
    
    for i in 0:N, j in 0:N
        push!(all_index_pairs, (i, j))
        push!(all_values, √(i^2 + j^2))
    end
    
    sorted_indices = sortperm(all_values)
    
    values = []
    index_pairs = []

    for i in 1:(N+1)^2
        if all_values[sorted_indices[i]] > N break end

        if all_values[sorted_indices[i]] in values
            push!(index_pairs[end], all_index_pairs[sorted_indices[i]])
        else
            push!(values, all_values[sorted_indices[i]])
            push!(index_pairs, [all_index_pairs[sorted_indices[i]]])
        end
    end

    return values, index_pairs
end


