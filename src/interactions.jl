"""
    count_interactions(;
        pdb_file,
        traj_file,
        fs_file,
        output_file,
        lipids,
        donors=[(lipid, lipid.head_atom) for lipid in lipids],
        acceptors=[(lipid, "PO4") for lipid in lipids],
        cutoff=6.0)

Counting the number of interactions (such as H-bonds) defined as a donor lipid atom being closer than the cutoff to a acceptor lipid atom. Saves the results in an HDF5 file with the label `tag`.

### Keyword arguments

* `pdb_file`: PDB structure file;
* `traj_files`: a list of trajectory files;
* `fs_files`: a list of corresponding HDF5 fluctuation spectrum files of the trajectory files;
* `output_file`: output directory;
* `lipids`: a list of lipids of type `Lipid` as defined in `lipids.jl`;
* `donors`: list of tuples of lipids (from type `Lipid`) and atoms to be considered as interaction donors;
* `acceptors`: list of tuples of lipids (from type `Lipid`) and atoms to be considered as interaction acceptors;
* `cutoff`: interaction distance cutoff.

"""
function count_interactions(;
    pdb_file,
    traj_file,
    fs_file,
    output_file,
    lipids,
    donors=[(lipid, lipid.head_atom) for lipid in lipids],
    acceptors=[(lipid, "PO4") for lipid in lipids],
    cutoff=6.0
)
    # finding indices of tail atoms

    atoms = PDBTools.readPDB(pdb_file)

    tail_atoms_inds = [a.index for a in atoms if (a.resname, a.name) in [(lipid.name, lipid.tail_atom) for lipid in lipids]]

    # determining donor and acceptor atoms indices

    donor_inds = []

    for donor in donors
        append!(donor_inds, [a.index for a in atoms if a.resname == donor[1].name && a.name == donor[2]])
    end

    acceptor_inds = []

    for acceptor in acceptors
        append!(acceptor_inds, [a.index for a in atoms if a.resname == acceptor[1].name && a.name == acceptor[2]])
    end

    # create index to resnum dictionary

    index_to_resnum = Dict()

    for a in atoms
        if a.resname in [lipid.name for lipid in lipids]
            index_to_resnum[a.index] = a.resnum
        end
    end

    # reading trajectory    

    traj = Chemfiles.Trajectory(traj_file)
    n_frames = Int(Chemfiles.size(traj))
    
    l_id = h5read(fs_file, "l_id")
    
    # create matrix to save number of interactions of each lipid in each frame
    
    n_lipids = length(unique([a.resnum for a in atoms if a.resname in [lipid.name for lipid in lipids]]))
    tag = zeros(Int, n_frames, n_lipids)
    
    println("Identifying interactions for trajectory file $(traj_file)")
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
        
        # counting the interactions of each donor
        
        for ind_A in donor_inds
            n_interactions = 0
            for ind_B in acceptor_inds
                if leaflet_id[ind_A] == leaflet_id[ind_B] && peuclidean(coords[:, ind_A], coords[:, ind_B], (Lx, Ly, Lz)) < cutoff && index_to_resnum[ind_A] != index_to_resnum[ind_B]
                    n_interactions += 1
                end
            end
            tag[frame_index, index_to_resnum[ind_A]] = n_interactions
        end
    end

    Chemfiles.close(traj); GC.gc()
    
    h5write(output_file, "tag", tag)
    
    return nothing
end

