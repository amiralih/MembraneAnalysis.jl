"""
    curvature(;
        point,
        hq,
        box_dims,
        q_max)

Calculates curvature per mode of a point from fluctuation spectrum.

### Keyword arguments

* `point`: ordered pair of X and Y values;
* `hq`: 2D matrix of height fluctuation spectrum;
* `box_dims`: ordered pair of simulation box Lx and Ly values;
* `q_max`: maximum q mode magnitude value to be used.

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
    lipids_sampled_curvature(;
        pdb_file,
        traj_files,
        fs_files,
        output_dir,
        lipids,
        q_max)

Calculates mean sampled curvature of heavy atoms of each lipid. Results for lipid "XXXX" will be stored in `XXXX_cs.dat` in the output directory.

### Keyword arguments

* `pdb_file`: PDB structure file;
* `traj_files`: a list of trajectory files;
* `fs_files`: a list of corresponding HDF5 fluctuation spectrum files of the trajectory files;
* `output_dir`: output directory;
* `lipids`: a list of lipids of type `Lipid` as defined in `lipids.jl`;
* `q_max`: maximum q mode magnitude value to be used.

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
        l_id = h5read(fs_file, "l_id")

        println("Calculating mean curvatures for trajectory file $(traj_file)")
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
    peptide_sampled_curvature(;
        pdb_file,
        traj_files,
        fs_files,
        output_dir,
        lipids,
        q_max)

Calculates mean sampled curvature of CA atoms of peptide residues. 

### Keyword arguments

* `pdb_file`: PDB structure file;
* `traj_files`: a list of trajectory files;
* `fs_files`: a list of corresponding HDF5 fluctuation spectrum files of the trajectory files;
* `output_file`: output file;
* `lipids`: a list of lipids of type `Lipid` as defined in `lipids.jl`;
* `n_residues`: number of residues in the peptide;
* `q_max`: maximum q mode magnitude value to be used.

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
        l_id = h5read(fs_file, "l_id")

        println("Calculating mean curvature for trajectory file $(traj_file)")
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
    lipids_curvature_spectrum(;
        pdb_file,
        traj_files,
        fs_files,
        output_dir,
        lipids,
        ref_atoms=Dict(lipid => lipid.ref_atom for lipid in lipids),
        q_max)

Calculates curvature spectrum of the lipids using their reference atom position. Assumes square (Lx = Ly) bilayer. Results for lipid "XXXX" will be stored in `XXXX_cqs.dat` in the output directory.

### Keyword arguments

* `pdb_file`: PDB structure file;
* `traj_files`: a list of trajectory files;
* `fs_files`: a list of corresponding HDF5 fluctuation spectrum files of the trajectory files;
* `output_dir`: output directory;
* `lipids`: a list of lipids of type `Lipid` as defined in `lipids.jl`;
* `ref_atoms`: a dictionary of reference atoms for each lipid;
* `q_max`: maximum q mode magnitude value to be used.

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
        l_id = h5read(fs_file, "l_id")

        println("Calculating curvature spectrum for trajectory file $(traj_file)")
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
    peptide_curvature_spectrum(;
        pdb_file,
        traj_files,
        fs_files,
        output_file,
        lipids,
        ref_residue),
        q_max)

Calculates curvature spectrum of the peptide using the CA atom of its reference residue. Assumes square (Lx = Ly) bilayer.

### Keyword arguments

* `pdb_file`: PDB structure file;
* `traj_files`: a list of trajectory files;
* `fs_files`: a list of corresponding HDF5 fluctuation spectrum files of the trajectory files;
* `output_file`: output file;
* `lipids`: a list of lipids of type `Lipid` as defined in `lipids.jl`;
* `ref_residue`: residue number of the reference residue of the peptide;
* `q_max`: maximum q mode magnitude value to be used.

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
        l_id = h5read(fs_file, "l_id")

        println("Calculating curvature spectrum for trajectory file $(traj_file)")
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

"""
    TCB_analysis(;
        input_dir,
        lipids,
        weights=ones(length(lipids)) ./ length(lipids),
        z_cutoff,
        area=readdlm(input_dir * "A.dat")[1],
        output_dir,
        tcb_plot=false)

Calculates bilayer bending rigidity modulus and mean sampled curvature of lipids relative to a weighted average from transverse curvature bias analysis. Optionally plots mean sampled curvature of atoms of each lipid as a function of height.

### Keyword arguments

* `input_dir`: directory with lipid atoms height and curvature files (e.g. `XXXX_zs.dat` and `XXXX_cs.dat` for lipid "XXXX");
* `lipids`: a list of lipids of type `Lipid` as defined in `lipids.jl`;
* `weights`: a list of the same size as `lipids`, determining the weight of each lipid's TCB curve in the analysis. (Should be equal to the fraction of bilayer area covered by that lipid. Will be equal by default.);
* `z_cutoff`: cutoff height to exclude anomalous behavior near lipid head region;
* `area`: bilayer area, will be read from `A.dat` in `input_dir` by default;
* `output_dir`: output directory;
* `tcb_plot`: saves a plot ("TCB_plot.pdf") in `output_dir` if `true`.

"""
function TCB_analysis(;
    input_dir,
    lipids,
    weights=ones(length(lipids)) ./ length(lipids),
    z_cutoff,
    area=readdlm(input_dir * "A.dat")[1],
    output_dir,
    tcb_plot=false
)
    atoms_z = Dict()
    atoms_c = Dict()
    z = Dict()
    c = Dict()
    c_e = Dict()
    
    @. model(x, p) = p[1] + p[2] * x
    p = Dict()

    # plot  <c>(z) vs z
    if tcb_plot
        default(fontfamily="Computer Modern", framestyle=:box,
                legend=:outerright, grid=false,
                xlabel=L"z",
                ylabel=L"\langle c \rangle(z-\delta)",
                yformatter=:scientific)
        plot()

        hline!([0], lc=:black, label=false)
    end

    slope = 0.0
    intercept = 0.0
    
    for (lipid, w) in zip(lipids, weights)
        data = readdlm(input_dir * "$(lipid.name)_zs.dat")
        z[lipid] = data[:, 2]
        z_order = sortperm(z[lipid])

        atoms_z[lipid] = data[:, 1]

        data = readdlm(input_dir * "$(lipid.name)_cs.dat")
        c[lipid] = data[:, 2]
        c_e[lipid] = data[:, 3]

        atoms_c[lipid] = data[:, 1]

        @assert atoms_z[lipid] == atoms_c[lipid]

        z[lipid] = z[lipid][z_order]
        c[lipid] = c[lipid][z_order]
        c_e[lipid] = c_e[lipid][z_order]

        zs = z[lipid]
        cs = c[lipid]
        fit = curve_fit(model, zs[zs .< z_cutoff], cs[zs .< z_cutoff], [0.0, 0.0])
        p[lipid] = coef(fit)

        slope += w * p[lipid][2]
        intercept += w * p[lipid][1]

        if tcb_plot
            scatter!(z[lipid], c[lipid], yerror=c_e[lipid],
                     msc=:auto, label=lipid.name)
        end
    end

    if tcb_plot
    vline!([z_cutoff], lc=:black, ls=:dash, label=false)
    savefig(output_dir * "TCB_plot.pdf")
    end
    
    kc = -2 / (area * slope)

    writedlm(output_dir * "kc.dat", [kc])
    writedlm(output_dir * "delta_cs.dat",
             [[lipid.name for lipid in lipids] [(p[lipid][1] - intercept) for lipid in lipids]])

    return nothing
end

