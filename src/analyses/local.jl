# Heights from the midplane, sampled curvature and curvature/thickness spectra of lipids and
# proteins.

function _lipid_species(system::MembraneSystem, lipids)
    lipids === nothing && return collect(eachindex(system.lipids))
    lipids isa AbstractString && (lipids = [lipids])
    isempty(lipids) && throw(ArgumentError("the list of lipids is empty"))
    species = [lipid_index(system, String(n)) for n in lipids]
    allunique(species) || throw(ArgumentError("lipids are listed more than once: $(collect(lipids))"))
    return species
end

# suffix of the analysis id of analyses restricted to some of the lipid species
_species_tag(system::MembraneSystem, species::Vector{Int}) =
    species == collect(eachindex(system.lipids)) ? "" :
    "." * join((system.lipids[s].name for s in species), "+")

function _protein_residue_list(system::MembraneSystem, p::Int, residues)
    residues === nothing && return protein_residues(system, p)
    isempty(residues) && throw(ArgumentError("the list of residues is empty"))
    return collect(Int, residues)
end

_lipid_files(a, system) = ["$(system.lipids[s].name)_$(a.label).dat" for s in a.species]

_direction_text(d::Symbol) = d === :all ? "all directions" : d === :x ? "modes along x" : "modes along y"

_check_direction(d::Symbol) =
    d in (:all, :x, :y) || throw(ArgumentError("direction must be :all, :x or :y (got :$d)"))

# ----------------------------------------------------------------------- lipid atom heights

"""
    LipidHeights(system; lipids=nothing, label="zs")

Mean height above the bilayer midplane of every height atom (see [`Lipid`](@ref)) of each lipid
species, signed so that it is positive towards the lipid's own leaflet surface.

Output: `<LIPID>_<label>.dat` with columns atom, mean [Å], sem.
"""
struct LipidHeights <: AbstractAnalysis
    label::String
    tag::String
    species::Vector{Int}
    groups::Vector{Vector{AtomGroup}}
    labels::Vector{String}
end

function LipidHeights(system::MembraneSystem; lipids=nothing, label::AbstractString="zs")
    species = _lipid_species(system, lipids)
    groups = [height_atom_groups(system, s) for s in species]
    labels = [system.lipids[s].name * ":" * g.label for (s, gs) in zip(species, groups) for g in gs]
    return LipidHeights(String(label), _species_tag(system, species), species, groups, labels)
end

analysis_id(a::LipidHeights) = "lipid_heights." * a.label * a.tag
output_files(a::LipidHeights, system) = _lipid_files(a, system)
component_labels(a::LipidHeights) = a.labels

function accumulate_frame!(sums, counts, a::LipidHeights, system, ctx::FrameContext)
    c = 0
    for gs in a.groups, g in gs
        c += 1
        _accumulate_heights!(sums, counts, c, g, system, ctx)
    end
    return nothing
end

function write_results(a::LipidHeights, bs::BlockStats, system, out::OutputContext)
    m = block_means(bs)
    e = block_sems(bs)
    c = 0
    for (s, gs) in zip(a.species, a.groups)
        r = (c + 1):(c + length(gs))
        c += length(gs)
        name = system.lipids[s].name
        write_table(out, "$(name)_$(a.label).dat",
                    "$name: mean height of atoms above the bilayer midplane [Å]; " * stats_header(bs),
                    ["atom", "mean", "sem"], [g.label for g in gs], m[r], e[r])
    end
    return Dict{String,Any}()
end

# ---------------------------------------------------------------------- protein heights

"""
    ProteinHeights(system; protein, residues=nothing, label="zs")

Mean height of the backbone atom(s) of every residue of a protein above the bilayer midplane
(sign convention of the protein's `leaflet`, see [`Protein`](@ref)). Copies are averaged.

Output: `<PROTEIN>_<label>.dat` with columns residue, mean [Å], sem.
"""
struct ProteinHeights <: AbstractAnalysis
    label::String
    protein::Int
    pname::String
    groups::Vector{AtomGroup}
    labels::Vector{String}
end

function ProteinHeights(system::MembraneSystem; protein::AbstractString, residues=nothing,
                        label::AbstractString="zs")
    p = protein_index(system, protein)
    groups = [protein_residue_group(system, p, r) for r in _protein_residue_list(system, p, residues)]
    labels = ["$(protein):$(g.label)" for g in groups]
    return ProteinHeights(String(label), p, String(protein), groups, labels)
end

analysis_id(a::ProteinHeights) = "protein_heights.$(a.pname)_$(a.label)"
output_files(a::ProteinHeights, system) = ["$(a.pname)_$(a.label).dat"]
component_labels(a::ProteinHeights) = a.labels

function accumulate_frame!(sums, counts, a::ProteinHeights, system, ctx::FrameContext)
    for (c, g) in enumerate(a.groups)
        _accumulate_heights!(sums, counts, c, g, system, ctx)
    end
    return nothing
end

function write_results(a::ProteinHeights, bs::BlockStats, system, out::OutputContext)
    name = system.proteins[a.protein].name
    write_table(out, "$(name)_$(a.label).dat",
                "$name: mean height of residue backbone atoms above the bilayer midplane [Å]; " *
                stats_header(bs), ["residue", "mean", "sem"], [g.label for g in a.groups],
                block_means(bs), block_sems(bs))
    return Dict{String,Any}()
end

# ------------------------------------------------------------------ lipid sampled curvature

"""
    LipidSampledCurvature(system; q_max, q_min=0.0, direction=:all, lipids=nothing, label="cs")

Mean curvature per mode sampled by every height atom of each lipid species, from the height
spectrum restricted to `q_min < |q| < q_max` (and to modes along x or y for
`direction = :x`/`:y`). Signed so that positive curvature bends towards the lipid's leaflet.

Output: `<LIPID>_<label>.dat` with columns atom, mean [1/Å], sem (`label` defaults to `cs`,
`cs_X` or `cs_Y`).
"""
struct LipidSampledCurvature <: AbstractAnalysis
    label::String
    tag::String
    species::Vector{Int}
    groups::Vector{Vector{AtomGroup}}
    shell::QBins
    direction::Symbol
    labels::Vector{String}
end

function LipidSampledCurvature(system::MembraneSystem; q_max::Real, q_min::Real=0.0,
                               direction::Symbol=:all, lipids=nothing,
                               label::AbstractString="cs" * direction_suffix(direction))
    _check_direction(direction)
    species = _lipid_species(system, lipids)
    groups = [height_atom_groups(system, s) for s in species]
    labels = [system.lipids[s].name * ":" * g.label for (s, gs) in zip(species, groups) for g in gs]
    return LipidSampledCurvature(String(label), _species_tag(system, species), species, groups,
                                 q_shell(q_min, q_max), direction, labels)
end

analysis_id(a::LipidSampledCurvature) = "lipid_sampled_curvature." * a.label * a.tag
output_files(a::LipidSampledCurvature, system) = _lipid_files(a, system)
component_labels(a::LipidSampledCurvature) = a.labels

function accumulate_frame!(sums, counts, a::LipidSampledCurvature, system, ctx::FrameContext)
    field = FourierField(ctx.hq, ctx.Lx, ctx.Ly, a.shell, a.direction, :curvature)
    c = 0
    for gs in a.groups, g in gs
        _accumulate_field!(sums, counts, c, field, g, system, ctx, true)
        c += 1
    end
    return nothing
end

function write_results(a::LipidSampledCurvature, bs::BlockStats, system, out::OutputContext)
    m = block_means(bs)
    e = block_sems(bs)
    c = 0
    for (s, gs) in zip(a.species, a.groups)
        r = (c + 1):(c + length(gs))
        c += length(gs)
        name = system.lipids[s].name
        write_table(out, "$(name)_$(a.label).dat",
                    "$name: mean sampled curvature per mode [1/Å], $(a.shell.lo[1]) < q < " *
                    "$(a.shell.hi[1]) 1/Å, $(_direction_text(a.direction)); " * stats_header(bs),
                    ["atom", "mean", "sem"], [g.label for g in gs], m[r], e[r])
    end
    return Dict{String,Any}()
end

# ---------------------------------------------------------------- protein sampled curvature

"""
    ProteinSampledCurvature(system; protein, q_max, q_min=0.0, direction=:all,
                            residues=nothing, label="cs")

Mean curvature per mode sampled by the backbone atom(s) of every residue of a protein.

Output: `<PROTEIN>_<label>.dat` with columns residue, mean [1/Å], sem.
"""
struct ProteinSampledCurvature <: AbstractAnalysis
    label::String
    protein::Int
    pname::String
    groups::Vector{AtomGroup}
    shell::QBins
    direction::Symbol
    labels::Vector{String}
end

function ProteinSampledCurvature(system::MembraneSystem; protein::AbstractString, q_max::Real,
                                 q_min::Real=0.0, direction::Symbol=:all, residues=nothing,
                                 label::AbstractString="cs" * direction_suffix(direction))
    _check_direction(direction)
    p = protein_index(system, protein)
    groups = [protein_residue_group(system, p, r) for r in _protein_residue_list(system, p, residues)]
    labels = ["$(protein):$(g.label)" for g in groups]
    return ProteinSampledCurvature(String(label), p, String(protein), groups,
                                   q_shell(q_min, q_max), direction, labels)
end

analysis_id(a::ProteinSampledCurvature) = "protein_sampled_curvature.$(a.pname)_$(a.label)"
output_files(a::ProteinSampledCurvature, system) = ["$(a.pname)_$(a.label).dat"]
component_labels(a::ProteinSampledCurvature) = a.labels

function accumulate_frame!(sums, counts, a::ProteinSampledCurvature, system, ctx::FrameContext)
    field = FourierField(ctx.hq, ctx.Lx, ctx.Ly, a.shell, a.direction, :curvature)
    for (c, g) in enumerate(a.groups)
        _accumulate_field!(sums, counts, c - 1, field, g, system, ctx, true)
    end
    return nothing
end

function write_results(a::ProteinSampledCurvature, bs::BlockStats, system, out::OutputContext)
    name = system.proteins[a.protein].name
    write_table(out, "$(name)_$(a.label).dat",
                "$name: mean sampled curvature per mode of residue backbone atoms [1/Å], " *
                "$(a.shell.lo[1]) < q < $(a.shell.hi[1]) 1/Å, $(_direction_text(a.direction)); " *
                stats_header(bs), ["residue", "mean", "sem"], [g.label for g in a.groups],
                block_means(bs), block_sems(bs))
    return Dict{String,Any}()
end

# ------------------------------------------------------------------------- lipid spectra

"""
    LipidSpectrum(system, kind; q_min, q_max, n_bins, direction=:all, lipids=nothing,
                  ref_atoms=Dict(), label)

Curvature (`kind = :curvature`) or thickness (`kind = :thickness`) spectrum of each lipid
species: the mean curvature (thickness deviation) per mode at the lipid's reference atom in
`n_bins` shells between `q_min` and `q_max`. Curvature is signed towards the lipid's leaflet.
`ref_atoms` maps lipid names to reference atoms that replace the lipid's `ref_atom`.

Output: `<LIPID>_<label>.dat` with columns q (bin centre) [1/Å], mean, sem. Bins without modes
give `NaN`. `label` defaults to `cqs`/`tqs` plus `_X`/`_Y` for directional spectra.
"""
struct LipidSpectrum <: AbstractAnalysis
    kind::Symbol
    label::String
    tag::String
    species::Vector{Int}
    groups::Vector{AtomGroup}
    bins::QBins
    direction::Symbol
    labels::Vector{String}
end

function LipidSpectrum(system::MembraneSystem, kind::Symbol; q_min::Real, q_max::Real,
                       n_bins::Integer, direction::Symbol=:all, lipids=nothing,
                       ref_atoms::AbstractDict=Dict{String,String}(),
                       label::AbstractString=(kind === :curvature ? "cqs" : "tqs") * direction_suffix(direction))
    kind in (:curvature, :thickness) || throw(ArgumentError("kind must be :curvature or :thickness"))
    _check_direction(direction)
    species = _lipid_species(system, lipids)
    for (name, atom) in ref_atoms
        lipid_index(system, String(name))   # error for unknown lipids
        atom isa AbstractString || throw(ArgumentError("ref_atoms: the atom of $name must be a string"))
    end
    groups = AtomGroup[]
    for s in species
        lipid = system.lipids[s]
        atom = String(get(ref_atoms, lipid.name, lipid.ref_atom))
        push!(groups, AtomGroup(atom, find_lipid_atoms(system.atoms, lipid.name, atom;
                                                       pdb_file=system.pdb_file), Int8[]))
    end
    bins = spectrum_bins(q_min, q_max, n_bins)
    labels = ["$(system.lipids[s].name):q$(b)" for s in species for b in 1:n_bins]
    return LipidSpectrum(kind, String(label), _species_tag(system, species), species, groups,
                         bins, direction, labels)
end

analysis_id(a::LipidSpectrum) = "lipid_$(a.kind)_spectrum." * a.label * a.tag
output_files(a::LipidSpectrum, system) = _lipid_files(a, system)
component_labels(a::LipidSpectrum) = a.labels

function accumulate_frame!(sums, counts, a::LipidSpectrum, system, ctx::FrameContext)
    spectrum = a.kind === :curvature ? ctx.hq : ctx.tq
    field = FourierField(spectrum, ctx.Lx, ctx.Ly, a.bins, a.direction, a.kind)
    nb = length(a.bins)
    for (k, g) in enumerate(a.groups)
        _accumulate_field!(sums, counts, (k - 1) * nb, field, g, system, ctx, a.kind === :curvature)
    end
    return nothing
end

function write_results(a::LipidSpectrum, bs::BlockStats, system, out::OutputContext)
    m = block_means(bs)
    e = block_sems(bs)
    nb = length(a.bins)
    what = a.kind === :curvature ? "curvature spectrum: mean curvature per mode [1/Å]" :
           "thickness spectrum: mean thickness deviation per mode [Å]"
    for (k, s) in enumerate(a.species)
        r = ((k - 1) * nb + 1):(k * nb)
        name = system.lipids[s].name
        write_table(out, "$(name)_$(a.label).dat",
                    "$name $what at reference atom $(a.groups[k].label), $(_direction_text(a.direction)); " *
                    stats_header(bs), ["q [1/Å]", "mean", "sem"], a.bins.centers, m[r], e[r])
    end
    return Dict{String,Any}()
end

# ----------------------------------------------------------------------- protein spectra

"""
    ProteinSpectrum(system, kind; protein, residue, q_min, q_max, n_bins, direction=:all, label)

Curvature or thickness spectrum of a protein at the backbone atom(s) of residue `residue`
(averaged over copies of the protein).

Output: `<PROTEIN>_<label>.dat` with columns q (bin centre) [1/Å], mean, sem. `label` defaults
to `cqs`/`tqs` plus `_X`/`_Y` for directional spectra.
"""
struct ProteinSpectrum <: AbstractAnalysis
    kind::Symbol
    label::String
    protein::Int
    pname::String
    group::AtomGroup
    bins::QBins
    direction::Symbol
    labels::Vector{String}
end

function ProteinSpectrum(system::MembraneSystem, kind::Symbol; protein::AbstractString,
                         residue::Integer, q_min::Real, q_max::Real, n_bins::Integer,
                         direction::Symbol=:all,
                         label::AbstractString=(kind === :curvature ? "cqs" : "tqs") * direction_suffix(direction))
    kind in (:curvature, :thickness) || throw(ArgumentError("kind must be :curvature or :thickness"))
    _check_direction(direction)
    p = protein_index(system, protein)
    group = protein_residue_group(system, p, residue)
    bins = spectrum_bins(q_min, q_max, n_bins)
    labels = ["$(protein):q$(b)" for b in 1:n_bins]
    return ProteinSpectrum(kind, String(label), p, String(protein), group, bins, direction, labels)
end

analysis_id(a::ProteinSpectrum) = "protein_$(a.kind)_spectrum.$(a.pname)_$(a.label)"
output_files(a::ProteinSpectrum, system) = ["$(a.pname)_$(a.label).dat"]
component_labels(a::ProteinSpectrum) = a.labels

function accumulate_frame!(sums, counts, a::ProteinSpectrum, system, ctx::FrameContext)
    spectrum = a.kind === :curvature ? ctx.hq : ctx.tq
    field = FourierField(spectrum, ctx.Lx, ctx.Ly, a.bins, a.direction, a.kind)
    _accumulate_field!(sums, counts, 0, field, a.group, system, ctx, a.kind === :curvature)
    return nothing
end

function write_results(a::ProteinSpectrum, bs::BlockStats, system, out::OutputContext)
    name = system.proteins[a.protein].name
    what = a.kind === :curvature ? "curvature spectrum: mean curvature per mode [1/Å]" :
           "thickness spectrum: mean thickness deviation per mode [Å]"
    write_table(out, "$(name)_$(a.label).dat",
                "$name $what at the backbone of residue $(a.group.label), " *
                "$(_direction_text(a.direction)); " * stats_header(bs),
                ["q [1/Å]", "mean", "sem"], a.bins.centers, block_means(bs), block_sems(bs))
    return Dict{String,Any}()
end
