"""
    fluctuation_spectrum(;
        pdb_file,
        traj_file,
        output_file,
        lipids,
        ref_atoms=Dict(lipid => lipid.ref_atom for lipid in lipids),
        L_grid)

Calculates the height and thickness fluctuation spectrum of a lipid bilayer simulation trajectory and saves the results as a HDF5 file with labels `hq` and `tq`.

### Keyword arguments

* `pdf_file`: PDB structure file;
* `traj_file`: trajectory file;
* `output_file`: output HDF5;
* `lipids`: a list of lipids of type `Lipid` as defined in `lipids.jl`;
* `ref_atoms`: a dictionary of reference atoms for each lipid;
* `L_grid`: length of the lattice grid used to discretize the surface.

"""
function fluctuation_spectrum(; 
    pdb_file,
    traj_file,
    output_file,
    lipids,
    ref_atoms=Dict(lipid => lipid.ref_atom for lipid in lipids),
    L_grid
)

    # finding indices of reference atoms

    atoms = PDBTools.readPDB(pdb_file)
    leaflet_id = atom_leaflet(pdb_file=pdb_file, lipids=lipids)

    tail_atoms_inds = [a.index for a in atoms if (a.resname, a.name) in [(lipid.name, lipid.tail_atom) for lipid in lipids]]

    ref_inds = []

    for lipid in lipids
        append!(ref_inds, [a.index for a in atoms if a.resname == lipid.name && a.name == ref_atoms[lipid]])
    end

    # calculating number of grid cells

    (Lx, Ly, _) = box_dimensions(traj_file=traj_file) 
    n_grid_x = Int(ceil(Lx / L_grid))
    n_grid_y = Int(ceil(Ly / L_grid))

    # opening trajectory

    traj = Chemfiles.Trajectory(traj_file)
    n_frames = Int(Chemfiles.size(traj))

    hq = zeros(ComplexF32, n_frames, n_grid_x, n_grid_y)
    tq = zeros(ComplexF32, n_frames, n_grid_x, n_grid_y)
    l_id = zeros(Int32, n_frames, length(leaflet_id))

    println("Calculating fluctuation spectrum for trajectory file $(traj_file)")
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
        
        ref_inds_1 = [i for i in ref_inds if leaflet_id[i] == 1]
        ref_inds_2 = [i for i in ref_inds if leaflet_id[i] == -1]
        
        # calculate initial mean Z values to dynamically determine lipids leaflet
        
        zs_m = zeros(Float32, n_grid_x, n_grid_y)

        # find average z of cells for leaflet 1

        zs_1 = zeros(Float32, n_grid_x, n_grid_y)
        ns_1 = zeros(Int32, n_grid_x, n_grid_y)
        
        for index in ref_inds_1
            x_index = Int(floor(coords[1, index] * n_grid_x / Lx)) + 1
            y_index = Int(floor(coords[2, index] * n_grid_y / Ly)) + 1   

            zs_1[x_index, y_index] += coords[3, index]
            ns_1[x_index, y_index] += 1
        end

        # assign weighted average of adjacent cells to empty cells

        for i in 1:n_grid_x, j in 1:n_grid_y
            if ns_1[i, j] == 0
                
                # find index pairs of adjacent cells

                up    = (mod1(i - 1, n_grid_x), j)
                down  = (mod1(i + 1, n_grid_x), j)
                left  = (i, mod1(j - 1, n_grid_y))
                right = (i, mod1(j + 1, n_grid_y))
                
                neighbors = [up, down, left, right]

                zs_1[i, j] = sum([zs_1[index_pair...] for index_pair in neighbors])
                ns_1[i, j] = sum([ns_1[index_pair...] for index_pair in neighbors])
            end
        end

        zs_1 = zs_1 ./ ns_1
        zs_m .+= zs_1
       
        # find average z of cells for leaflet 2

        zs_2 = zeros(Float32, n_grid_x, n_grid_y)
        ns_2 = zeros(Int32, n_grid_x, n_grid_y)
        
        for index in ref_inds_2
            x_index = Int(floor(coords[1, index] * n_grid_x / Lx)) + 1
            y_index = Int(floor(coords[2, index] * n_grid_y / Ly)) + 1   

            zs_2[x_index, y_index] += coords[3, index]
            ns_2[x_index, y_index] += 1
        end

        # assign weighted average of adjacent cells to empty cells

        for i in 1:n_grid_x, j in 1:n_grid_y
            if ns_2[i, j] == 0
                
                # find index pairs of adjacent cells

                up    = (mod1(i - 1, n_grid_x), j)
                down  = (mod1(i + 1, n_grid_x), j)
                left  = (i, mod1(j - 1, n_grid_y))
                right = (i, mod1(j + 1, n_grid_y))
                
                neighbors = [up, down, left, right]

                zs_2[i, j] = sum([zs_2[index_pair...] for index_pair in neighbors])
                ns_2[i, j] = sum([ns_2[index_pair...] for index_pair in neighbors])
            end
        end

        zs_2 = zs_2 ./ ns_2
        zs_m .+= zs_2

        # determine lipids leaflet
        
        zs_m ./= 2
        leaflet_id = atom_leaflet_dynamic(coords=coords, zs_m=zs_m, Ls=(Lx, Ly),
                                  pdb_file=pdb_file, lipids=lipids)
        l_id[frame_index, :] = leaflet_id
       
        # repeat analysis with dynamically determined lipid leaflet
        
        ref_inds_1 = [i for i in ref_inds if leaflet_id[i] == 1]
        ref_inds_2 = [i for i in ref_inds if leaflet_id[i] == -1]

        # find average z of cells for leaflet 1

        zs_1 = zeros(Float32, n_grid_x, n_grid_y)
        ns_1 = zeros(Int32, n_grid_x, n_grid_y)
        
        for index in ref_inds_1
            x_index = Int(floor(coords[1, index] * n_grid_x / Lx)) + 1
            y_index = Int(floor(coords[2, index] * n_grid_y / Ly)) + 1   

            zs_1[x_index, y_index] += coords[3, index]
            ns_1[x_index, y_index] += 1
        end

        # assign weighted average of adjacent cells to empty cells

        for i in 1:n_grid_x, j in 1:n_grid_y
            if ns_1[i, j] == 0
                
                # find index pairs of adjacent cells

                up    = (mod1(i - 1, n_grid_x), j)
                down  = (mod1(i + 1, n_grid_x), j)
                left  = (i, mod1(j - 1, n_grid_y))
                right = (i, mod1(j + 1, n_grid_y))
                
                neighbors = [up, down, left, right]

                zs_1[i, j] = sum([zs_1[index_pair...] for index_pair in neighbors])
                ns_1[i, j] = sum([ns_1[index_pair...] for index_pair in neighbors])
            end
        end

        zs_1 = zs_1 ./ ns_1
        zs_1 = zs_1 .- midplane
       
        # find average z of cells for leaflet 2

        zs_2 = zeros(Float32, n_grid_x, n_grid_y)
        ns_2 = zeros(Int32, n_grid_x, n_grid_y)
        
        for index in ref_inds_2
            x_index = Int(floor(coords[1, index] * n_grid_x / Lx)) + 1
            y_index = Int(floor(coords[2, index] * n_grid_y / Ly)) + 1   

            zs_2[x_index, y_index] += coords[3, index]
            ns_2[x_index, y_index] += 1
        end

        # assign weighted average of adjacent cells to empty cells

        for i in 1:n_grid_x, j in 1:n_grid_y
            if ns_2[i, j] == 0
                
                # find index pairs of adjacent cells

                up    = (mod1(i - 1, n_grid_x), j)
                down  = (mod1(i + 1, n_grid_x), j)
                left  = (i, mod1(j - 1, n_grid_y))
                right = (i, mod1(j + 1, n_grid_y))
                
                neighbors = [up, down, left, right]

                zs_2[i, j] = sum([zs_2[index_pair...] for index_pair in neighbors])
                ns_2[i, j] = sum([ns_2[index_pair...] for index_pair in neighbors])
            end
        end

        zs_2 = zs_2 ./ ns_2
        zs_2 = zs_2 .- midplane
        
        # calculate leaflets sepectra

        hq_1 = fft(zs_1) * ((Lx * Ly) / (n_grid_x * n_grid_y))
        hq_2 = fft(zs_2) * ((Lx * Ly) / (n_grid_x * n_grid_y))

        hq[frame_index, :, :]  = (hq_1 .+ hq_2) ./ 2
        tq[frame_index, :, :]  = (hq_1 .- hq_2)
        
        # phase-shifting the FFT to correspond to cell centers

        Δϕs_x = fftfreq(n_grid_x, n_grid_x/Lx) * (-π) * (Lx / n_grid_x)
        Δϕs_y = fftfreq(n_grid_y, n_grid_y/Ly) * (-π) * (Ly / n_grid_y)
        Δϕs = [(Δϕ_x + Δϕ_y) for Δϕ_x in Δϕs_x, Δϕ_y in Δϕs_y]

        hq[frame_index, :, :] .*= exp.(im .* Δϕs)
        tq[frame_index, :, :] .*= exp.(im .* Δϕs)
    end
    
    Chemfiles.close(traj); GC.gc()
    
    h5write(output_file, "hq", hq)
    h5write(output_file, "tq", tq)
    h5write(output_file, "l_id", l_id)

    return nothing
end

"""
    area_expansion_modulus(;
        traj_files,
        output_file)

Calculates area expansion modulus in units of kBT per square of length unit in trajectories (e.g., Å^2).

### Keyword arguments

* `traj_files`: a list of trajectory files;
* `output_file`: output file to save the result.

"""
function area_expansion_modulus(;
        traj_files,
        output_file
)
    As = Float64[]

    for traj_file in traj_files

        traj = Chemfiles.Trajectory(traj_file)
        n_frames = Int(Chemfiles.size(traj))

        println("Finding simulation box area for trajectory file $(traj_file)")
        @showprogress 0.01 "Analyzing trajectory..." for frame_index in 1:n_frames
            frame = Chemfiles.read_step(traj, frame_index - 1)
            box_dims = Chemfiles.lengths(Chemfiles.UnitCell(frame))
            push!(As, box_dims[1] * box_dims[2])
        end

        Chemfiles.close(traj); GC.gc()
    end
    
    A_m = mean(As)
    A_e = blocking_error(As)
    
    VarA_m = mean((As .- A_m).^2)
    VarA_e = blocking_error((As .- A_m).^2)
    
    KA_m = A_m / VarA_m    
    KA_e = KA_m * √((A_e / A_m)^2 + (VarA_e / VarA_m)^2)
    
    writedlm(output_file, [KA_m KA_e])

    return (KA_m, KA_e)
end

"""
    fluctuation_hq2_data(;
        box_dims,
        fs_files,
        output_dir,
        q_max)

Calculates |hq|^2 vs |q| data from fluctuation spectrum. Assumes square (Lx = Ly) bilayer. Saves results to `hq2.dat` in `output_dir`.

### Keyword arguments

* `box_dims`: ordered pair of lateral simulation box dimensions, Lx and Ly;
* `fs_files`: a list of HDF5 fluctuation spectrum files;
* `output_dir`: output directory;
* `q_max`: maximum q mode magnitude value to be used.

"""
function fluctuation_hq2_data(;
	box_dims,
	fs_files,
	output_dir,
    q_max
)
	# find |q| values (up to q_max) and corresponing indices

	(L, _) = box_dims

	# including at least 9 q values (for N = 4) even if greater than q_max

	N = max(Int(floor(q_max / (2π / L))), 4)
	values, index_pairs = get_index_pairs(N)
	qs = values * (2π / L)

	# make lists of |hq|^2 values (for each |q|) of all frames

	hq2s = zeros(length(qs) - 1)
	hq2es = zeros(length(qs) - 1)

	for q_i in 2:length(qs)
		hq2 = Float64[]
		for ind_pair in index_pairs[q_i]
	    		for fs_file in fs_files
				append!(hq2, (abs.(h5read(fs_file, "hq")[:,  (ind_pair .+ 1)...])).^2)
			end
		end
		
		hq2s[q_i - 1] = mean(hq2)
		hq2es[q_i - 1] = blocking_error(hq2)
	end

	writedlm(output_dir * "hq2.dat", [qs[2:end] hq2s hq2es])

	return nothing
end
