"""
Calculates the lateral density spectrum of lipids in a lipid bilayer simulation trajectory and saves the output as a HDF5 file.
"""
function lipids_density_spectrum(; 
    pdb_file,
    traj_file,
    fs_file,
    output_file,
    lipids,
    ref_atoms=Dict(lipid => lipid.ref_atom for lipid in lipids),
    L_grid
)

    # finding indices of reference atoms

    atoms = PDBTools.readPDB(pdb_file)
    leaflet_id = atom_leaflet(pdb_file=pdb_file, lipids=lipids)
    l_id = h5read(fs_file, "l_id")

    tail_atoms_inds = [a.index for a in atoms if (a.resname, a.name) in [(lipid.name, lipid.tail_atom) for lipid in lipids]]

    ref_inds = Dict()
    ref_inds_all = []

    for lipid in lipids
        ref_inds[lipid.name] = [a.index for a in atoms if a.resname == lipid.name && a.name == ref_atoms[lipid]]
        append!(ref_inds_all, [a.index for a in atoms if a.resname == lipid.name && a.name == ref_atoms[lipid]])
    end

    # calculating number of grid cells

    (Lx, Ly, _) = box_dimensions(traj_file=traj_file) 
    n_grid_x = Int(ceil(Lx / L_grid))
    n_grid_y = Int(ceil(Ly / L_grid))

    # opening trajectory

    traj = Chemfiles.Trajectory(traj_file)
    n_frames = Int(Chemfiles.size(traj))

    dq_1 = Dict()
    dq_2 = Dict()

    for lipid in lipids
        dq_1[lipid.name] = zeros(ComplexF32, n_frames, n_grid_x, n_grid_y)
        dq_2[lipid.name] = zeros(ComplexF32, n_frames, n_grid_x, n_grid_y)
    end

    println("Calculating density spectrum for trajectory file $(traj_file)")
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
        
        for lipid in lipids
                
            ref_inds_1 = [i for i in ref_inds[lipid.name] if leaflet_id[i] == 1]
            ref_inds_2 = [i for i in ref_inds[lipid.name] if leaflet_id[i] == -1]
            
            # find the number in each cell for leaflets

            ns_1 = zeros(Int32, n_grid_x, n_grid_y)
            
            for index in ref_inds_1
                x_index = Int(floor(coords[1, index] * n_grid_x / Lx)) + 1
                y_index = Int(floor(coords[2, index] * n_grid_y / Ly)) + 1   

                ns_1[x_index, y_index] += 1
            end

            ns_2 = zeros(Int32, n_grid_x, n_grid_y)
            
            for index in ref_inds_2
                x_index = Int(floor(coords[1, index] * n_grid_x / Lx)) + 1
                y_index = Int(floor(coords[2, index] * n_grid_y / Ly)) + 1   

                ns_2[x_index, y_index] += 1
            end

            # calculate leaflets sepectra

            dq_1[lipid.name][frame_index, :, :] = fft(ns_1) * ((Lx * Ly) / (n_grid_x * n_grid_y))
            dq_2[lipid.name][frame_index, :, :] = fft(ns_2) * ((Lx * Ly) / (n_grid_x * n_grid_y))

            # phase-shifting the FFT to correspond to cell centers

            Δϕs_x = fftfreq(n_grid_x, n_grid_x/Lx) * (-π) * (Lx / n_grid_x)
            Δϕs_y = fftfreq(n_grid_y, n_grid_y/Ly) * (-π) * (Ly / n_grid_y)
            Δϕs = [(Δϕ_x + Δϕ_y) for Δϕ_x in Δϕs_x, Δϕ_y in Δϕs_y]

            dq_1[lipid.name][frame_index, :, :] .*= exp.(im .* Δϕs)
            dq_2[lipid.name][frame_index, :, :] .*= exp.(im .* Δϕs)
        end

    end
    
    Chemfiles.close(traj); GC.gc()
   
    for lipid in lipids
        h5write(output_file, lipid.name * "_1", dq_1[lipid.name])
        h5write(output_file, lipid.name * "_2", dq_2[lipid.name])
    end

    return nothing
end

"""
Calculates the lateral density spectrum of the peptide in a lipid bilayer simulation trajectory and saves the output as a HDF5 file.
"""
function peptide_density_spectrum(; 
    pdb_file,
    traj_file,
    fs_file,
    output_file,
    lipids,
    ref_residue,
    L_grid
)

    # finding indices of reference atoms of lipids

    atoms = PDBTools.readPDB(pdb_file)
    leaflet_id = atom_leaflet(pdb_file=pdb_file, lipids=lipids)
    l_id = h5read(fs_file, "l_id")

    tail_atoms_inds = [a.index for a in atoms if (a.resname, a.name) in [(lipid.name, lipid.tail_atom) for lipid in lipids]]

    # finding reference atoms of the peptide in each leaflet

    ref_inds_pep = [a.index for a in atoms if a.segname[1:3] == "PRO" && a.resnum == ref_residue && a.name == "CA"]
    ref_inds_1 = [j for j in ref_inds_pep if atoms[j].z > 0.0]
    ref_inds_2 = [j for j in ref_inds_pep if atoms[j].z < 0.0]
    
    # calculating number of grid cells

    (Lx, Ly, _) = box_dimensions(traj_file=traj_file) 
    n_grid_x = Int(ceil(Lx / L_grid))
    n_grid_y = Int(ceil(Ly / L_grid))

    # opening trajectory

    traj = Chemfiles.Trajectory(traj_file)
    n_frames = Int(Chemfiles.size(traj))

    dq_1 = zeros(ComplexF32, n_frames, n_grid_x, n_grid_y)
    dq_2 = zeros(ComplexF32, n_frames, n_grid_x, n_grid_y)

    println("Calculating density spectrum for trajectory file $(traj_file)")
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
        
        # find the number in each cell for leaflets

        ns_1 = zeros(Int32, n_grid_x, n_grid_y)
        
        for index in ref_inds_1
            x_index = Int(floor(coords[1, index] * n_grid_x / Lx)) + 1
            y_index = Int(floor(coords[2, index] * n_grid_y / Ly)) + 1   

            ns_1[x_index, y_index] += 1
        end

        ns_2 = zeros(Int32, n_grid_x, n_grid_y)
        
        for index in ref_inds_2
            x_index = Int(floor(coords[1, index] * n_grid_x / Lx)) + 1
            y_index = Int(floor(coords[2, index] * n_grid_y / Ly)) + 1   

            ns_2[x_index, y_index] += 1
        end

        # calculate leaflets sepectra

        dq_1[frame_index, :, :] = fft(ns_1) * ((Lx * Ly) / (n_grid_x * n_grid_y))
        dq_2[frame_index, :, :] = fft(ns_2) * ((Lx * Ly) / (n_grid_x * n_grid_y))
        
        # phase-shifting the FFT to correspond to cell centers

        Δϕs_x = fftfreq(n_grid_x, n_grid_x/Lx) * (-π) * (Lx / n_grid_x)
        Δϕs_y = fftfreq(n_grid_y, n_grid_y/Ly) * (-π) * (Ly / n_grid_y)
        Δϕs = [(Δϕ_x + Δϕ_y) for Δϕ_x in Δϕs_x, Δϕ_y in Δϕs_y]

        dq_1[frame_index, :, :] .*= exp.(im .* Δϕs)
        dq_2[frame_index, :, :] .*= exp.(im .* Δϕs)
    end
    
    Chemfiles.close(traj); GC.gc()
   
    h5write(output_file, "PEP_1", dq_1)
    h5write(output_file, "PEP_2", dq_2)

    return nothing
end

"""
Calculates the lateral (2-D) radial density function of lipid A and lipid B in a lipid bilayer simulation trajectory.
"""
function lipids_radial_distribution(; 
    pdb_file,
    traj_file,
    fs_file,
    output_file,
    lipids,
    ref_atoms=Dict(lipid => lipid.ref_atom for lipid in lipids),
    lipid_A,
    lipid_B,
    max_r=20,
    n_bins=200,
    same_leaflet=true
)

    # finding indices of reference atoms

    atoms = PDBTools.readPDB(pdb_file)
    leaflet_id = atom_leaflet(pdb_file=pdb_file, lipids=lipids)
    l_id = h5read(fs_file, "l_id")

    tail_atoms_inds = [a.index for a in atoms if (a.resname, a.name) in [(lipid.name, lipid.tail_atom) for lipid in lipids]]

    ref_inds_A1 = [a.index for a in atoms if a.resname == lipid_A.name
                   && a.name == ref_atoms[lipid_A] && leaflet_id[a.index] == 1]
    ref_inds_A2 = [a.index for a in atoms if a.resname == lipid_A.name
                   && a.name == ref_atoms[lipid_A] && leaflet_id[a.index] == -1]
    ref_inds_B1 = [a.index for a in atoms if a.resname == lipid_B.name
                   && a.name == ref_atoms[lipid_B] && leaflet_id[a.index] == 1]
    ref_inds_B2 = [a.index for a in atoms if a.resname == lipid_B.name
                   && a.name == ref_atoms[lipid_B] && leaflet_id[a.index] == -1]

    # number of distances in bins

    bins = zeros(Int64, n_bins)

    # opening trajectory

    traj = Chemfiles.Trajectory(traj_file)
    n_frames = Int(Chemfiles.size(traj))

    # record area

    area = zeros(Float64, n_frames)

    println("Calculating radial density function for trajectory file $(traj_file)")
    @showprogress 0.01 "Analyzing trajectory..." for frame_index in 1:n_frames
        
        # read a frame

        leaflet_id = l_id[frame_index, :]

        ref_inds_A1 = [a.index for a in atoms if a.resname == lipid_A.name
                       && a.name == ref_atoms[lipid_A] && leaflet_id[a.index] == 1]
        ref_inds_A2 = [a.index for a in atoms if a.resname == lipid_A.name
                       && a.name == ref_atoms[lipid_A] && leaflet_id[a.index] == -1]
        ref_inds_B1 = [a.index for a in atoms if a.resname == lipid_B.name
                       && a.name == ref_atoms[lipid_B] && leaflet_id[a.index] == 1]
        ref_inds_B2 = [a.index for a in atoms if a.resname == lipid_B.name
                       && a.name == ref_atoms[lipid_B] && leaflet_id[a.index] == -1]

        frame = Chemfiles.read_step(traj, frame_index - 1)
        box_dims = Chemfiles.lengths(Chemfiles.UnitCell(frame))
        Lx, Ly, Lz = box_dims
        coords = Chemfiles.positions(frame)

        area[frame_index] = Lx * Ly

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
        
        # calculating distances for all selected atoms

        if same_leaflet
            for atom_1 in ref_inds_A1
            for atom_2 in ref_inds_B1
                d = peuclidean(coords[1:2, atom_1], coords[1:2, atom_2], (Lx, Ly))
                if atom_1 != atom_2 && d < max_r
                    bin_index = Int(floor((d / max_r) * n_bins)) + 1
                    bins[bin_index] += 1
                end
            end
            end

            for atom_1 in ref_inds_A2
            for atom_2 in ref_inds_B2
                d = peuclidean(coords[1:2, atom_1], coords[1:2, atom_2], (Lx, Ly))
                if atom_1 != atom_2 && d < max_r
                    bin_index = Int(floor((d / max_r) * n_bins)) + 1
                    bins[bin_index] += 1
                end
            end
            end
        else
            for atom_1 in ref_inds_A1
            for atom_2 in ref_inds_B2
                d = peuclidean(coords[1:2, atom_1], coords[1:2, atom_2], (Lx, Ly))
                if atom_1 != atom_2 && d < max_r
                    bin_index = Int(floor((d / max_r) * n_bins)) + 1
                    bins[bin_index] += 1
                end
            end
            end

            for atom_1 in ref_inds_A2
            for atom_2 in ref_inds_B1
                d = peuclidean(coords[1:2, atom_1], coords[1:2, atom_2], (Lx, Ly))
                if atom_1 != atom_2 && d < max_r
                    bin_index = Int(floor((d / max_r) * n_bins)) + 1
                    bins[bin_index] += 1
                end
            end
            end
        end
    end
    
    Chemfiles.close(traj); GC.gc()
   
    bin_width = max_r / n_bins
    bin_centers = collect(range(bin_width/2, step=bin_width, length=n_bins))

    area = mean(area)
    
    N_A1 = length(ref_inds_A1)
    N_A2 = length(ref_inds_A2)
    N_B1 = length(ref_inds_B1)
    N_B2 = length(ref_inds_B2)

    if same_leaflet
        if lipid_A == lipid_B
            norm_factor = n_frames * ((N_A1 * (N_B1 - 1)) + (N_A2 * (N_B2 - 1))) * 2π * bin_width / area
        else
            norm_factor = n_frames * ((N_A1 * N_B1) + (N_A2 * N_B2)) * 2π * bin_width / area
        end
    else
        norm_factor = n_frames * ((N_A1 * N_B2) + (N_A2 * N_B1)) * 2π * bin_width / area
    end

    bins = bins ./ (norm_factor .* bin_centers)
    
    writedlm(output_file, [bin_centers bins])
    
    return nothing
end

"""
Calculates the lateral (2-D) radial density function of peptide and a lipid in a lipid bilayer simulation trajectory.
"""
function peptide_radial_distribution(; 
    pdb_file,
    traj_file,
    output_file,
    lipids,
    ref_atoms=Dict(lipid => lipid.ref_atom for lipid in lipids),
    ref_residue,
    lipid,
    max_r=20,
    n_bins=200,
    same_leaflet=true
)

    # finding indices of reference atoms

    atoms = PDBTools.readPDB(pdb_file)
    leaflet_id = atom_leaflet(pdb_file=pdb_file, lipids=lipids)
    l_id = h5read(fs_file, "l_id")

    tail_atoms_inds = [a.index for a in atoms if (a.resname, a.name) in [(lipid.name, lipid.tail_atom) for lipid in lipids]]
    
    ref_inds_A1 = [a.index for a in atoms if a.segname[1:3] == "PRO"
                   && a.resnum == ref_residue && a.name == "CA" && a.z > 0]
    ref_inds_A2 = [a.index for a in atoms if a.segname[1:3] == "PRO"
                   && a.resnum == ref_residue && a.name == "CA" && a.z < 0]
    
    ref_inds_B1 = [a.index for a in atoms if a.resname == lipid.name
                   && a.name == ref_atoms[lipid] && leaflet_id[a.index] == 1]
    ref_inds_B2 = [a.index for a in atoms if a.resname == lipid.name
                   && a.name == ref_atoms[lipid] && leaflet_id[a.index] == -1]

    # number of distances in bins

    bins = zeros(Int64, n_bins)

    # opening trajectory

    traj = Chemfiles.Trajectory(traj_file)
    n_frames = Int(Chemfiles.size(traj))

    # record area

    area = zeros(Float64, n_frames)

    println("Calculating radial density function for trajectory file $(traj_file)")
    @showprogress 0.01 "Analyzing trajectory..." for frame_index in 1:n_frames
        
        # read a frame
        
        leaflet_id = l_id[frame_index, :]
    
        ref_inds_B1 = [a.index for a in atoms if a.resname == lipid.name
                       && a.name == ref_atoms[lipid] && leaflet_id[a.index] == 1]
        ref_inds_B2 = [a.index for a in atoms if a.resname == lipid.name
                       && a.name == ref_atoms[lipid] && leaflet_id[a.index] == -1]

        frame = Chemfiles.read_step(traj, frame_index - 1)
        box_dims = Chemfiles.lengths(Chemfiles.UnitCell(frame))
        Lx, Ly, Lz = box_dims
        coords = Chemfiles.positions(frame)

        area[frame_index] = Lx * Ly

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
        
        # calculating distances for all selected atoms

        if same_leaflet
            for atom_1 in ref_inds_A1
            for atom_2 in ref_inds_B1
                d = peuclidean(coords[1:2, atom_1], coords[1:2, atom_2], (Lx, Ly))
                if atom_1 != atom_2 && d < max_r
                    bin_index = Int(floor((d / max_r) * n_bins)) + 1
                    bins[bin_index] += 1
                end
            end
            end

            for atom_1 in ref_inds_A2
            for atom_2 in ref_inds_B2
                d = peuclidean(coords[1:2, atom_1], coords[1:2, atom_2], (Lx, Ly))
                if atom_1 != atom_2 && d < max_r
                    bin_index = Int(floor((d / max_r) * n_bins)) + 1
                    bins[bin_index] += 1
                end
            end
            end
        else
            for atom_1 in ref_inds_A1
            for atom_2 in ref_inds_B2
                d = peuclidean(coords[1:2, atom_1], coords[1:2, atom_2], (Lx, Ly))
                if atom_1 != atom_2 && d < max_r
                    bin_index = Int(floor((d / max_r) * n_bins)) + 1
                    bins[bin_index] += 1
                end
            end
            end

            for atom_1 in ref_inds_A2
            for atom_2 in ref_inds_B1
                d = peuclidean(coords[1:2, atom_1], coords[1:2, atom_2], (Lx, Ly))
                if atom_1 != atom_2 && d < max_r
                    bin_index = Int(floor((d / max_r) * n_bins)) + 1
                    bins[bin_index] += 1
                end
            end
            end
        end
    end
    
    Chemfiles.close(traj); GC.gc()
   
    bin_width = max_r / n_bins
    bin_centers = collect(range(bin_width/2, step=bin_width, length=n_bins))

    area = mean(area)
    
    N_A1 = length(ref_inds_A1)
    N_A2 = length(ref_inds_A2)
    N_B1 = length(ref_inds_B1)
    N_B2 = length(ref_inds_B2)

    if same_leaflet
        norm_factor = n_frames * ((N_A1 * N_B1) + (N_A2 * N_B2)) * 2π * bin_width / area
    else
        norm_factor = n_frames * ((N_A1 * N_B2) + (N_A2 * N_B1)) * 2π * bin_width / area
    end

    bins = bins ./ (norm_factor .* bin_centers)
    
    writedlm(output_file, [bin_centers bins])
    
    return nothing
end

