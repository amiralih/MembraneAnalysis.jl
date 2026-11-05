# The analysed system: atoms from the PDB file plus the lipid/protein bookkeeping that every
# analysis needs. Built once, before any trajectory is read.

"""
    AtomGroup(label, atoms, signs)

A set of atoms analysed together (e.g. all `C1A` atoms of POPC, or the backbone atoms of
residue 4 of a peptide). `signs` holds the fixed leaflet sign (±1) of every atom for proteins;
for lipids it is empty because their leaflet is determined anew in every frame.
"""
struct AtomGroup
    label::String
    atoms::Vector{Int}
    signs::Vector{Int8}
end

"""
    MembraneSystem(pdb_file, lipids, proteins=Protein[])

Read `pdb_file` and set up everything needed to analyse trajectories of this system:
lipid molecules and their leaflets in the structure, tail atoms (bilayer midplane), reference
atoms (membrane surface) and protein atoms. Inconsistencies between the lipid/protein
definitions and the structure (missing atoms, wrong `n_atoms`, …) raise an error.
"""
struct MembraneSystem
    pdb_file::String
    atoms::PDBAtoms
    n_atoms::Int
    lipids::Vector{Lipid}
    proteins::Vector{Protein}
    tail_atoms::Vector{Int}           # lipid tail atoms; their mean height is the midplane
    mol_head::Vector{Int}             # first atom of every lipid molecule
    mol_species::Vector{Int}          # lipid species (index into `lipids`) of every molecule
    atom_mol::Vector{Int32}           # lipid molecule of every atom (0: not a lipid atom)
    static_leaflet::Vector{Int8}      # leaflet (±1) of every molecule in the PDB structure
    ref_atoms::Vector{Vector{Int}}    # reference atoms of every lipid species
    surface_atoms::Vector{Int}        # all reference atoms (species order), define the surface
    pdb_midplane::Float64             # midplane height in the PDB structure
    protein_atoms::Vector{Vector{Int}}     # all atoms of every protein
    protein_backbone::Vector{Vector{Int}}  # backbone atoms of every protein
    protein_sign::Vector{Int8}             # leaflet sign of every protein atom (0 otherwise)
end

function MembraneSystem(pdb_file::AbstractString, lipids::AbstractVector{Lipid},
                        proteins::AbstractVector{Protein}=Protein[])
    isempty(lipids) && throw(ArgumentError("at least one lipid species must be defined"))
    names = [l.name for l in lipids]
    allunique(names) || throw(ArgumentError("lipid names must be unique, got $(names)"))
    pnames = [p.name for p in proteins]
    allunique(pnames) || throw(ArgumentError("protein names must be unique, got $(pnames)"))

    atoms = read_pdb(pdb_file)
    N = length(atoms)

    # lipid molecules and their leaflet in the structure
    mol_head = Int[]; mol_species = Int[]; static_leaflet = Int8[]
    atom_mol = zeros(Int32, N)
    for (s, lipid) in enumerate(lipids)
        heads = [i for i in 1:N if atoms.resname[i] == lipid.name && atoms.name[i] == lipid.head_atom]
        isempty(heads) && throw(ArgumentError(
            "lipid $(lipid.name): no atom named \"$(lipid.head_atom)\" in residues named " *
            "\"$(lipid.name)\" in $pdb_file"))
        for h in heads
            last = h + lipid.n_atoms - 1
            if last > N || atoms.residue[last] != atoms.residue[h] ||
               (last < N && atoms.residue[last + 1] == atoms.residue[h])
                throw(ArgumentError(
                    "lipid $(lipid.name): n_atoms = $(lipid.n_atoms) does not match the molecule " *
                    "starting at atom $h (residue $(atoms.resnum[h])) in $pdb_file. The head atom " *
                    "must be the first atom of the lipid and n_atoms its number of atoms."))
            end
            if h > 1 && atoms.residue[h - 1] == atoms.residue[h]
                throw(ArgumentError(
                    "lipid $(lipid.name): head atom \"$(lipid.head_atom)\" is not the first atom " *
                    "of residue $(atoms.resnum[h]) in $pdb_file"))
            end
            push!(mol_head, h)
            push!(mol_species, s)
            atom_mol[h:last] .= length(mol_head)
            push!(static_leaflet, atoms.z[h] > atoms.z[last] ? Int8(1) : Int8(-1))
        end
    end

    lipid_names = Set(names)
    stray = findfirst(i -> atom_mol[i] == 0 && atoms.resname[i] in lipid_names, 1:N)
    if stray !== nothing
        throw(ArgumentError(
            "atom $stray ($(atoms.resname[stray]) $(atoms.name[stray]), residue " *
            "$(atoms.resnum[stray])) in $pdb_file is not part of any lipid molecule; check the " *
            "head_atom and n_atoms of lipid $(atoms.resname[stray])"))
    end

    pairs = Set((l.name, l.tail_atom) for l in lipids)
    tail_atoms = [i for i in 1:N if (atoms.resname[i], atoms.name[i]) in pairs]
    for l in lipids
        any(i -> atoms.resname[i] == l.name, tail_atoms) || throw(ArgumentError(
            "lipid $(l.name): no tail atom named \"$(l.tail_atom)\" found in $pdb_file"))
    end

    ref_atoms = [find_lipid_atoms(atoms, l.name, l.ref_atom; pdb_file=pdb_file) for l in lipids]
    surface_atoms = reduce(vcat, ref_atoms)
    pdb_midplane = mean(atoms.z[tail_atoms])

    # proteins
    protein_atoms = Vector{Int}[]
    protein_backbone = Vector{Int}[]
    protein_sign = zeros(Int8, N)
    for p in proteins
        sel = [i for i in 1:N if _in_protein(atoms, p, i)]
        isempty(sel) && throw(ArgumentError("protein $(p.name): no atoms match its selection " *
                                            "($(p)) in $pdb_file"))
        bb = [i for i in sel if atoms.name[i] in p.backbone_atoms]
        isempty(bb) && throw(ArgumentError("protein $(p.name): no backbone atoms named " *
                                           "$(p.backbone_atoms) found in $pdb_file"))
        for i in sel
            protein_sign[i] = p.leaflet === :upper ? Int8(1) :
                              p.leaflet === :lower ? Int8(-1) :
                              (atoms.z[i] > pdb_midplane ? Int8(1) : Int8(-1))
        end
        push!(protein_atoms, sel)
        push!(protein_backbone, bb)
    end

    return MembraneSystem(String(pdb_file), atoms, N, collect(lipids), collect(proteins),
                          tail_atoms, mol_head, mol_species, atom_mol, static_leaflet,
                          ref_atoms, surface_atoms, pdb_midplane, protein_atoms,
                          protein_backbone, protein_sign)
end

function _in_protein(atoms::PDBAtoms, p::Protein, i::Int)
    if !isempty(p.segments)
        return atoms.segname[i] in p.segments
    elseif !isempty(p.chains)
        return atoms.chain[i] in p.chains
    else
        return startswith(atoms.segname[i], p.segment_prefix)
    end
end

function find_lipid_atoms(atoms::PDBAtoms, resname::AbstractString, atomname::AbstractString;
                          pdb_file::AbstractString="the PDB file")
    inds = [i for i in 1:length(atoms) if atoms.resname[i] == resname && atoms.name[i] == atomname]
    isempty(inds) && throw(ArgumentError(
        "no atom named \"$atomname\" in residues named \"$resname\" in $pdb_file"))
    return inds
end

Base.show(io::IO, s::MembraneSystem) =
    print(io, "MembraneSystem(\"", basename(s.pdb_file), "\": ", s.n_atoms, " atoms, ",
          length(s.mol_head), " lipids of ", length(s.lipids), " species, ",
          length(s.proteins), " protein(s))")

function Base.show(io::IO, ::MIME"text/plain", s::MembraneSystem)
    println(io, "MembraneSystem from ", s.pdb_file)
    println(io, "  atoms: ", s.n_atoms)
    for (k, l) in enumerate(s.lipids)
        n = count(==(k), s.mol_species)
        n_up = count(m -> s.mol_species[m] == k && s.static_leaflet[m] == 1, eachindex(s.mol_head))
        println(io, "  ", l, ": ", n, " molecules (", n_up, " upper / ", n - n_up, " lower)")
    end
    for (k, p) in enumerate(s.proteins)
        res = protein_residues(s, k)
        println(io, "  ", p, ": ", length(s.protein_atoms[k]), " atoms, residues ",
                isempty(res) ? "-" : "$(first(res))–$(last(res))")
    end
    print(io, "  bilayer midplane in structure: z = ", round(s.pdb_midplane; digits=2), " Å")
end

# --------------------------------------------------------------------------- lookups

"""
    lipid_index(system, name) -> Int

Index of the lipid species called `name`.
"""
function lipid_index(s::MembraneSystem, name::AbstractString)
    k = findfirst(l -> l.name == name, s.lipids)
    k === nothing && throw(ArgumentError("unknown lipid \"$name\"; defined lipids: " *
                                         join((l.name for l in s.lipids), ", ")))
    return k
end

"""
    protein_index(system, name) -> Int

Index of the protein called `name`.
"""
function protein_index(s::MembraneSystem, name::AbstractString)
    k = findfirst(p -> p.name == name, s.proteins)
    k === nothing && throw(ArgumentError("unknown protein \"$name\"; defined proteins: " *
                                         (isempty(s.proteins) ? "none" :
                                          join((p.name for p in s.proteins), ", "))))
    return k
end

"""
    height_atom_groups(system, k) -> Vector{AtomGroup}

Atom groups (one per atom name, in PDB order) of lipid species `k` used by the height,
sampled-curvature and TCB analyses.
"""
function height_atom_groups(s::MembraneSystem, k::Int)
    lipid = s.lipids[k]
    a = s.atoms
    if isempty(lipid.height_atoms)
        names = String[]
        for i in 1:s.n_atoms
            if a.resname[i] == lipid.name && !isempty(a.name[i]) &&
               first(a.name[i]) in ('C', 'D', 'G', 'X') && !(a.name[i] in names)
                push!(names, a.name[i])
            end
        end
        isempty(names) && throw(ArgumentError(
            "lipid $(lipid.name): no atoms with names starting with C, D, G or X; " *
            "set `height_atoms` for this lipid"))
    else
        names = lipid.height_atoms
    end
    return [AtomGroup(n, find_lipid_atoms(a, lipid.name, n; pdb_file=s.pdb_file), Int8[])
            for n in names]
end

"""
    protein_residues(system, k) -> Vector{Int}

Residue numbers (in order of appearance) that have backbone atoms in protein `k`.
"""
protein_residues(s::MembraneSystem, k::Int) = unique(s.atoms.resnum[s.protein_backbone[k]])

"""
    protein_residue_group(system, k, residue) -> AtomGroup

Backbone atoms of residue number `residue` in all copies of protein `k`, with their fixed
leaflet signs.
"""
function protein_residue_group(s::MembraneSystem, k::Int, residue::Integer)
    inds = [i for i in s.protein_backbone[k] if s.atoms.resnum[i] == residue]
    isempty(inds) && throw(ArgumentError(
        "protein $(s.proteins[k].name): residue $residue has no backbone atoms " *
        "$(s.proteins[k].backbone_atoms)"))
    return AtomGroup(string(residue), inds, s.protein_sign[inds])
end
