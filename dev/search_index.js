var documenterSearchIndex = {"docs":
[{"location":"thickness_analysis.html#Thickness-Analysis","page":"Thickness Analysis","title":"Thickness Analysis","text":"","category":"section"},{"location":"thickness_analysis.html","page":"Thickness Analysis","title":"Thickness Analysis","text":"lipids_thickness_spectrum","category":"page"},{"location":"thickness_analysis.html#MembraneAnalysis.lipids_thickness_spectrum","page":"Thickness Analysis","title":"MembraneAnalysis.lipids_thickness_spectrum","text":"lipids_thickness_spectrum(;\n    pdb_file,\n    traj_files,\n    fs_files,\n    output_dir,\n    lipids,\n    ref_atoms=Dict(lipid => lipid.ref_atom for lipid in lipids),\n    q_max)\n\nCalculates thickness spectrum of the lipids using their reference atom position. Assumes square (Lx = Ly) bilayer. Results for lipid \"XXXX\" will be stored in XXXX_tqs.dat in the output directory.\n\nKeyword arguments\n\npdb_file: PDB structure file;\ntraj_files: a list of trajectory files;\nfs_files: a list of corresponding HDF5 fluctuation spectrum files of the trajectory files;\noutput_dir: output directory;\nlipids: a list of lipids of type Lipid as defined in lipids.jl;\nref_atoms: a dictionary of reference atoms for each lipid;\nq_max: maximum q mode magnitude value to be used.\n\n\n\n\n\n","category":"function"},{"location":"thickness_analysis.html","page":"Thickness Analysis","title":"Thickness Analysis","text":"peptide_thickness_spectrum","category":"page"},{"location":"thickness_analysis.html#MembraneAnalysis.peptide_thickness_spectrum","page":"Thickness Analysis","title":"MembraneAnalysis.peptide_thickness_spectrum","text":"peptide_thickness_spectrum(;\n    pdb_file,\n    traj_files,\n    fs_files,\n    output_file,\n    lipids,\n    ref_residue),\n    q_max)\n\nCalculates thickness spectrum of the peptide using the CA atom of its reference residue. Assumes square (Lx = Ly) bilayer.\n\nKeyword arguments\n\npdb_file: PDB structure file;\ntraj_files: a list of trajectory files;\nfs_files: a list of corresponding HDF5 fluctuation spectrum files of the trajectory files;\noutput_file: output file;\nlipids: a list of lipids of type Lipid as defined in lipids.jl;\nref_residue: residue number of the reference residue of the peptide;\nq_max: maximum q mode magnitude value to be used.\n\n\n\n\n\n","category":"function"},{"location":"curvature_analysis.html#Curvature-Analysis","page":"Curvature Analysis","title":"Curvature Analysis","text":"","category":"section"},{"location":"curvature_analysis.html","page":"Curvature Analysis","title":"Curvature Analysis","text":"curvature","category":"page"},{"location":"curvature_analysis.html#MembraneAnalysis.curvature","page":"Curvature Analysis","title":"MembraneAnalysis.curvature","text":"curvature(;\n    point,\n    hq,\n    box_dims,\n    q_max)\n\nCalculates curvature per mode of a point from fluctuation spectrum.\n\nKeyword arguments\n\npoint: ordered pair of X and Y values;\nhq: 2D matrix of height fluctuation spectrum;\nbox_dims: ordered pair of simulation box Lx and Ly values;\nq_max: maximum q mode magnitude value to be used.\n\n\n\n\n\n","category":"function"},{"location":"curvature_analysis.html","page":"Curvature Analysis","title":"Curvature Analysis","text":"lipids_sampled_curvature","category":"page"},{"location":"curvature_analysis.html#MembraneAnalysis.lipids_sampled_curvature","page":"Curvature Analysis","title":"MembraneAnalysis.lipids_sampled_curvature","text":"lipids_sampled_curvature(;\n    pdb_file,\n    traj_files,\n    fs_files,\n    output_dir,\n    lipids,\n    q_max)\n\nCalculates mean sampled curvature of heavy atoms of each lipid. Results for lipid \"XXXX\" will be stored in XXXX_cs.dat in the output directory.\n\nKeyword arguments\n\npdb_file: PDB structure file;\ntraj_files: a list of trajectory files;\nfs_files: a list of corresponding HDF5 fluctuation spectrum files of the trajectory files;\noutput_dir: output directory;\nlipids: a list of lipids of type Lipid as defined in lipids.jl;\nq_max: maximum q mode magnitude value to be used.\n\n\n\n\n\n","category":"function"},{"location":"curvature_analysis.html","page":"Curvature Analysis","title":"Curvature Analysis","text":"peptide_sampled_curvature","category":"page"},{"location":"curvature_analysis.html#MembraneAnalysis.peptide_sampled_curvature","page":"Curvature Analysis","title":"MembraneAnalysis.peptide_sampled_curvature","text":"peptide_sampled_curvature(;\n    pdb_file,\n    traj_files,\n    fs_files,\n    output_dir,\n    lipids,\n    q_max)\n\nCalculates mean sampled curvature of CA atoms of peptide residues. \n\nKeyword arguments\n\npdb_file: PDB structure file;\ntraj_files: a list of trajectory files;\nfs_files: a list of corresponding HDF5 fluctuation spectrum files of the trajectory files;\noutput_file: output file;\nlipids: a list of lipids of type Lipid as defined in lipids.jl;\nn_residues: number of residues in the peptide;\nq_max: maximum q mode magnitude value to be used.\n\n\n\n\n\n","category":"function"},{"location":"curvature_analysis.html","page":"Curvature Analysis","title":"Curvature Analysis","text":"lipids_curvature_spectrum","category":"page"},{"location":"curvature_analysis.html#MembraneAnalysis.lipids_curvature_spectrum","page":"Curvature Analysis","title":"MembraneAnalysis.lipids_curvature_spectrum","text":"lipids_curvature_spectrum(;\n    pdb_file,\n    traj_files,\n    fs_files,\n    output_dir,\n    lipids,\n    ref_atoms=Dict(lipid => lipid.ref_atom for lipid in lipids),\n    q_max)\n\nCalculates curvature spectrum of the lipids using their reference atom position. Assumes square (Lx = Ly) bilayer. Results for lipid \"XXXX\" will be stored in XXXX_cqs.dat in the output directory.\n\nKeyword arguments\n\npdb_file: PDB structure file;\ntraj_files: a list of trajectory files;\nfs_files: a list of corresponding HDF5 fluctuation spectrum files of the trajectory files;\noutput_dir: output directory;\nlipids: a list of lipids of type Lipid as defined in lipids.jl;\nref_atoms: a dictionary of reference atoms for each lipid;\nq_max: maximum q mode magnitude value to be used.\n\n\n\n\n\n","category":"function"},{"location":"curvature_analysis.html","page":"Curvature Analysis","title":"Curvature Analysis","text":"peptide_curvature_spectrum","category":"page"},{"location":"curvature_analysis.html#MembraneAnalysis.peptide_curvature_spectrum","page":"Curvature Analysis","title":"MembraneAnalysis.peptide_curvature_spectrum","text":"peptide_curvature_spectrum(;\n    pdb_file,\n    traj_files,\n    fs_files,\n    output_file,\n    lipids,\n    ref_residue),\n    q_max)\n\nCalculates curvature spectrum of the peptide using the CA atom of its reference residue. Assumes square (Lx = Ly) bilayer.\n\nKeyword arguments\n\npdb_file: PDB structure file;\ntraj_files: a list of trajectory files;\nfs_files: a list of corresponding HDF5 fluctuation spectrum files of the trajectory files;\noutput_file: output file;\nlipids: a list of lipids of type Lipid as defined in lipids.jl;\nref_residue: residue number of the reference residue of the peptide;\nq_max: maximum q mode magnitude value to be used.\n\n\n\n\n\n","category":"function"},{"location":"curvature_analysis.html","page":"Curvature Analysis","title":"Curvature Analysis","text":"TCB_analysis","category":"page"},{"location":"curvature_analysis.html#MembraneAnalysis.TCB_analysis","page":"Curvature Analysis","title":"MembraneAnalysis.TCB_analysis","text":"TCB_analysis(;\n    input_dir,\n    lipids,\n    weights=ones(length(lipids)) ./ length(lipids),\n    z_cutoff,\n    area=readdlm(input_dir * \"A.dat\")[1],\n    output_dir,\n    tcb_plot=false)\n\nCalculates bilayer bending rigidity modulus and mean sampled curvature of lipids relative to a weighted average from transverse curvature bias analysis. Optionally plots mean sampled curvature of atoms of each lipid as a function of height.\n\nKeyword arguments\n\ninput_dir: directory with lipid atoms height and curvature files (e.g. XXXX_zs.dat and XXXX_cs.dat for lipid \"XXXX\");\nlipids: a list of lipids of type Lipid as defined in lipids.jl;\nweights: a list of the same size as lipids, determining the weight of each lipid's TCB curve in the analysis. (Should be equal to the fraction of bilayer area covered by that lipid. Will be equal by default.);\nz_cutoff: cutoff height to exclude anomalous behavior near lipid head region;\narea: bilayer area, will be read from A.dat in input_dir by default;\noutput_dir: output directory;\ntcb_plot: saves a plot (TCB_plot.pdf) in output_dir if true.\n\n\n\n\n\n","category":"function"},{"location":"fluctuation_analysis.html#Fluctuation-Analysis","page":"Fluctuation Analysis","title":"Fluctuation Analysis","text":"","category":"section"},{"location":"fluctuation_analysis.html","page":"Fluctuation Analysis","title":"Fluctuation Analysis","text":"fluctuation_spectrum","category":"page"},{"location":"fluctuation_analysis.html#MembraneAnalysis.fluctuation_spectrum","page":"Fluctuation Analysis","title":"MembraneAnalysis.fluctuation_spectrum","text":"fluctuation_spectrum(;\n    pdb_file,\n    traj_file,\n    output_file,\n    lipids,\n    ref_atoms=Dict(lipid => lipid.ref_atom for lipid in lipids),\n    L_grid)\n\nCalculates the height and thickness fluctuation spectrum of a lipid bilayer simulation trajectory and saves the results as a HDF5 file with labels hq and tq.\n\nKeyword arguments\n\npdf_file: PDB structure file;\ntraj_file: trajectory file;\noutput_file: output HDF5;\nlipids: a list of lipids of type Lipid as defined in lipids.jl;\nref_atoms: a dictionary of reference atoms for each lipid;\nL_grid: length of the lattice grid used to discretize the surface.\n\n\n\n\n\n","category":"function"},{"location":"fluctuation_analysis.html","page":"Fluctuation Analysis","title":"Fluctuation Analysis","text":"area_expansion_modulus","category":"page"},{"location":"fluctuation_analysis.html#MembraneAnalysis.area_expansion_modulus","page":"Fluctuation Analysis","title":"MembraneAnalysis.area_expansion_modulus","text":"area_expansion_modulus(;\n    traj_files,\n    output_file)\n\nCalculates area expansion modulus in units of kBT per square of length unit in trajectories (e.g., Å^2).\n\nKeyword arguments\n\ntraj_files: a list of trajectory files;\noutput_file: output file to save the result.\n\n\n\n\n\n","category":"function"},{"location":"index.html#MembraneAnalysis.jl","page":"Introduction","title":"MembraneAnalysis.jl","text":"","category":"section"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"Documentation for MembraneAnalysis.jl package","category":"page"},{"location":"tutorial.html#Example-usage-of-MembraneAnalysis.jl","page":"Tutorial","title":"Example usage of MembraneAnalysis.jl","text":"","category":"section"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"In this tutorial we showcase some of the main functionality of the package by analyzing an atomistic simulations of a membrane containing DOPC lipids and %10 mol fraction cholesterol. MembraneAnalysis.jl needs a PDF structure file (here s.pdb) and trajectory files of any format supported by Chemfiles.jl (here t1.nc to t200.nc).","category":"page"},{"location":"tutorial.html#Setting-analysis-parameters","page":"Tutorial","title":"Setting analysis parameters","text":"","category":"section"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"It is useful to define a set of parameters that we will need for input arguments of the methods.","category":"page"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"pdb_file: the PDB structure file\ntraj_dir: folder containing trajectory files\ntraj_name: trajecotry file name without numerical index and file extension\ntraj_inds: range of numerical indices of trajectory files\ntraj_ext: file extension of trajecotries\nlipids: list of lipids in the system (from instances of Lipid type defined in src/lipids.jl\nL_grid: lattice grid length (in length unit of trajectories) used for discretizing the membrane surface\nq_max: q mode magnitude upper cut-off in the analysis\noutput_dir: a directory for all the output files generated by methods","category":"page"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"begin\n  const pdb_file = \"s.pdb\"\n  const traj_dir = \"./\"\n  const traj_name = \"t\"\n  const traj_inds = 1:200\n  const traj_ext = \".nc\"\n  const lipids = [DOPC_aa, Chol_aa]\n  const L_grid = 15\n  const q_max = 0.08\n  const output_dir = \"./\"\n  if !isdir(output_dir) mkdir(output_dir) end\nend","category":"page"},{"location":"tutorial.html#Calculating-fluctuation-spectrum","page":"Tutorial","title":"Calculating fluctuation spectrum","text":"","category":"section"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"fluctuation_spectrum calculates discrete Fourirer modes of surface height and thickness fluctuations and stores them in an HDF5 file for each trajectory.","category":"page"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"for i in traj_inds\n    traj_file = traj_dir * traj_name * string(i) * traj_ext\n    output_file = output_dir * \"fs_$(i).h5\"\n\n    fluctuation_spectrum(;\n        pdb_file=pdb_file,\n        traj_file=traj_file,\n        output_file=output_file,\n        lipids=lipids,\n        L_grid=L_grid\n    )\nend","category":"page"},{"location":"tutorial.html#Calculating-area-expansion-modulus","page":"Tutorial","title":"Calculating area expansion modulus","text":"","category":"section"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"area_expansion_modulus calculates area expansion modulus from area fluctuations and stores them in a file.","category":"page"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"traj_files = [traj_dir * traj_name * string(i) * traj_ext for i in traj_inds]\noutput_file = output_dir * \"KA.dat\"\n\narea_expansion_modulus(;\n    traj_files=traj_files,\n    output_file=output_file\n)","category":"page"},{"location":"tutorial.html#Calculating-mean-height-from-midplane-of-the-heavy-atoms-of-the-lipids","page":"Tutorial","title":"Calculating mean height from midplane of the heavy atoms of the lipids","text":"","category":"section"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"We use the first trajecotry to calculate the average height of the heavy atoms of the lipids using lipids_atoms_heights method which will save them to XXXX_zs.dat files in the specified output directory (where \"XXXX\" is the name of the lipid).","category":"page"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"traj_file = traj_dir * traj_name * string(first(traj_inds)) * traj_ext\nfs_file = output_dir * \"fs_$(first(traj_inds)).h5\"\n\nlipids_atoms_height(;\n    pdb_file=pdb_file,\n    traj_file=traj_file,\n    fs_file=fs_file,\n    output_dir=output_dir,\n    lipids=lipids\n)","category":"page"},{"location":"tutorial.html#Calculating-mean-sampled-curvature-of-heavy-atoms-of-the-lipids","page":"Tutorial","title":"Calculating mean sampled curvature of heavy atoms of the lipids","text":"","category":"section"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"We calculate the average sampled curvature of the heavy atoms of the lipids using lipids_sampled_curvature method which will save them to XXXX_cs.dat files in the specified output directory (where \"XXXX\" is the name of the lipid).","category":"page"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"traj_files = [traj_dir * traj_name * string(i) * traj_ext for i in traj_id]\nfs_files = [output_dir * \"fs_$(i).h5\" for i in traj_inds]\n\nlipids_sampled_curvature(;\n    pdb_file=pdb_file,\n    traj_files=traj_files,\n    fs_files=fs_files,\n    output_dir=output_dir,\n    lipids=lipids,\n    q_max=q_max\n)","category":"page"},{"location":"tutorial.html#Transverse-curvature-bias-curve","page":"Tutorial","title":"Transverse curvature bias curve","text":"","category":"section"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"Using the outputs of the two previous steps we can plot the transverse curvature bias curves, which can be used to identify the neutral surface atom of each lipid, calculate the bending modulus of the membrane and estimate the difference in spontanous curvature of the membrane components [1].","category":"page"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"<p align=\"center\">\n<img src=\"https://github.com/amiralih/MembraneAnalysis.jl/blob/64a2e086d41f8e4db37dca6d8748e8273798b6fa/docs/src/TCB_DOPC_10.png\" width=\"500\">\n</p>","category":"page"},{"location":"tutorial.html#References","page":"Tutorial","title":"References","text":"","category":"section"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"<a id=\"1\">[1]</a>  Sapp, K. C., Beaven, A. H., and Sodt, A. J. (2021).  Spatial extent of a single lipid's influence on bilayer mechanics Phys. Rev. E 103, 042413.","category":"page"}]
}
