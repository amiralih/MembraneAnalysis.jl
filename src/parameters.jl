# Parameters file (TOML) parsing and validation.

"""
    ParameterError(message)

Error in a parameters file or in the command-line arguments.
"""
struct ParameterError <: Exception
    msg::String
end
Base.showerror(io::IO, e::ParameterError) = print(io, "ParameterError: ", e.msg)

"""
    AnalysisConfig

Parsed (but not yet resolved against a structure) settings of a complete run. Create it with
[`load_parameters`](@ref) or [`parse_parameters`](@ref).
"""
struct AnalysisConfig
    title::String
    lipids::Vector{Lipid}
    proteins::Vector{Protein}
    grid_spacing::Float64
    grid_size::Union{Nothing,Tuple{Int,Int}}
    fill::Symbol
    mask_proteins::Vector{String}
    mask_residues::Union{Nothing,Vector{Int}}
    mask_radius::Float64
    center_threshold::Float64
    spectrum_defaults::Dict{String,Any}
    blocking::Symbol
    n_blocks::Int
    output_directory::String
    save_spectra::Bool
    analyses::Vector{Dict{String,Any}}
    source::String
end

const _TOP_KEYS = ["title", "lipids", "proteins", "surface", "spectrum", "statistics",
                   "output", "analysis"]
const _LIPID_KEYS = ["name", "preset", "head_atom", "ref_atom", "tail_atom", "n_atoms",
                     "height_atoms"]
const _PROTEIN_KEYS = ["name", "segments", "chains", "segment_prefix", "backbone_atoms",
                       "leaflet"]
const _SURFACE_KEYS = ["grid_spacing", "grid_size", "fill", "mask_proteins", "mask_residues",
                       "mask_radius", "center_threshold"]
const _SPECTRUM_KEYS = ["q_min", "q_max", "n_bins"]
const _STATISTICS_KEYS = ["blocking", "n_blocks"]
const _OUTPUT_KEYS = ["directory", "save_spectra"]

# keys accepted by every analysis type (besides "type" and "label")
const _ANALYSIS_KEYS = Dict(
    "box" => String[],
    "undulation" => ["q_max", "model", "n_fit"],
    "lipid_heights" => ["lipids"],
    "lipid_sampled_curvature" => ["lipids", "q_min", "q_max", "direction"],
    "lipid_curvature_spectrum" => ["lipids", "ref_atoms", "q_min", "q_max", "n_bins", "direction"],
    "lipid_thickness_spectrum" => ["lipids", "ref_atoms", "q_min", "q_max", "n_bins", "direction"],
    "protein_heights" => ["protein", "residues"],
    "protein_sampled_curvature" => ["protein", "residues", "q_min", "q_max", "direction"],
    "protein_curvature_spectrum" => ["protein", "residue", "q_min", "q_max", "n_bins", "direction"],
    "protein_thickness_spectrum" => ["protein", "residue", "q_min", "q_max", "n_bins", "direction"],
    "tcb" => ["lipids", "q_min", "q_max", "direction", "z_cutoff", "weights"],
)

"""
    analysis_types() -> Vector{String}

Names of the analysis types that can be requested in a parameters file.
"""
analysis_types() = sort!(collect(keys(_ANALYSIS_KEYS)))

function _check_keys(tbl::AbstractDict, allowed, loc::AbstractString)
    for k in keys(tbl)
        k in allowed || throw(ParameterError(
            "unknown key \"$k\" in $loc. Allowed keys: " * join(allowed, ", ")))
    end
    return nothing
end

function _get(tbl::AbstractDict, key::AbstractString, ::Type{T}, default, loc::AbstractString) where {T}
    haskey(tbl, key) || return default
    v = tbl[key]
    try
        if T === String
            v isa AbstractString || error("invalid value")
            return String(v)
        elseif T === Float64
            v isa Real || error("invalid value")
            return Float64(v)
        elseif T === Int
            (v isa Integer || (v isa Real && isinteger(v))) || error("invalid value")
            return Int(v)
        elseif T === Bool
            v isa Bool || error("invalid value")
            return v
        elseif T === Vector{String}
            v isa AbstractString && return [String(v)]
            v isa AbstractVector || error("invalid value")
            return String[String(x) for x in v]
        elseif T === Vector{Int}
            return _parse_int_list(v)
        else
            return convert(T, v)
        end
    catch
        throw(ParameterError("invalid value $(repr(v)) for \"$key\" in $loc (expected $(_type_text(T)))"))
    end
end

function _require(tbl::AbstractDict, key::AbstractString, ::Type{T}, loc::AbstractString) where {T}
    haskey(tbl, key) || throw(ParameterError("missing required key \"$key\" in $loc"))
    return _get(tbl, key, T, nothing, loc)
end

_type_text(::Type{String}) = "a string"
_type_text(::Type{Float64}) = "a number"
_type_text(::Type{Int}) = "an integer"
_type_text(::Type{Bool}) = "true or false"
_type_text(::Type{Vector{String}}) = "a list of strings"
_type_text(::Type{Vector{Int}}) = "a list of integers or a range such as \"1:8\""
_type_text(T) = string(T)

# [1, 2, 3], 4, "1:8", "1-8", "1,3,5"
function _parse_int_list(v)
    v isa Integer && return [Int(v)]
    if v isa AbstractVector
        all(x -> x isa Integer, v) || error("invalid value")
        return Int[x for x in v]
    end
    v isa AbstractString || error("invalid value")
    out = Int[]
    for part in split(v, ',')
        part = strip(part)
        m = match(r"^(-?\d+)\s*[:\-]\s*(-?\d+)$", part)
        if m !== nothing
            append!(out, parse(Int, m[1]):parse(Int, m[2]))
        else
            push!(out, parse(Int, part))
        end
    end
    return out
end

function _parse_direction(v, loc)
    v === nothing && return :all
    v isa AbstractString ||
        throw(ParameterError("invalid direction $(repr(v)) in $loc (use \"all\", \"x\" or \"y\")"))
    d = lowercase(String(v))
    d in ("all", "xy", "both") && return :all
    d == "x" && return :x
    d == "y" && return :y
    throw(ParameterError("invalid direction \"$v\" in $loc (use \"all\", \"x\" or \"y\")"))
end

function _parse_lipid(tbl::AbstractDict, k::Int)
    loc = "[[lipids]] entry $k"
    _check_keys(tbl, _LIPID_KEYS, loc)
    preset = _get(tbl, "preset", String, nothing, loc)
    if preset !== nothing
        haskey(LIPID_PRESETS, preset) || throw(ParameterError(
            "unknown lipid preset \"$preset\" in $loc. Available presets: " *
            join(sort!(collect(keys(LIPID_PRESETS))), ", ")))
        base = LIPID_PRESETS[preset]
        return Lipid(_get(tbl, "name", String, base.name, loc),
                     _get(tbl, "head_atom", String, base.head_atom, loc),
                     _get(tbl, "ref_atom", String, base.ref_atom, loc),
                     _get(tbl, "tail_atom", String, base.tail_atom, loc),
                     _get(tbl, "n_atoms", Int, base.n_atoms, loc),
                     _get(tbl, "height_atoms", Vector{String}, base.height_atoms, loc))
    end
    return Lipid(_require(tbl, "name", String, loc), _require(tbl, "head_atom", String, loc),
                 _require(tbl, "ref_atom", String, loc), _require(tbl, "tail_atom", String, loc),
                 _require(tbl, "n_atoms", Int, loc),
                 _get(tbl, "height_atoms", Vector{String}, String[], loc))
end

function _parse_protein(tbl::AbstractDict, k::Int)
    loc = "[[proteins]] entry $k"
    _check_keys(tbl, _PROTEIN_KEYS, loc)
    leaflet = _get(tbl, "leaflet", String, "auto", loc)
    leaflet in ("auto", "upper", "lower") ||
        throw(ParameterError("invalid leaflet \"$leaflet\" in $loc (use \"auto\", \"upper\" or \"lower\")"))
    return Protein(_require(tbl, "name", String, loc);
                   segments=_get(tbl, "segments", Vector{String}, String[], loc),
                   chains=_get(tbl, "chains", Vector{String}, String[], loc),
                   segment_prefix=_get(tbl, "segment_prefix", String, "PRO", loc),
                   backbone_atoms=_get(tbl, "backbone_atoms", Vector{String}, ["BB", "CA"], loc),
                   leaflet=Symbol(leaflet))
end

function _table(raw::AbstractDict, key::AbstractString)
    haskey(raw, key) || return Dict{String,Any}()
    v = raw[key]
    v isa AbstractDict || throw(ParameterError("\"$key\" must be a table ([$key])"))
    return v
end

function _table_array(raw::AbstractDict, key::AbstractString)
    haskey(raw, key) || return Dict{String,Any}[]
    v = raw[key]
    (v isa AbstractVector && all(x -> x isa AbstractDict, v)) ||
        throw(ParameterError("\"$key\" must be an array of tables ([[$key]])"))
    return v
end

"""
    parse_parameters(raw::AbstractDict; source="") -> AnalysisConfig

Validate a parsed TOML dictionary and fill in defaults. See the documentation for all
available keys.
"""
function parse_parameters(raw::AbstractDict; source::AbstractString="")
    _check_keys(raw, _TOP_KEYS, "the parameters file (top level)")
    title = _get(raw, "title", String, "", "the top level")

    lipid_tables = _table_array(raw, "lipids")
    isempty(lipid_tables) && throw(ParameterError("no lipids defined; add at least one [[lipids]] table"))
    lipids = [_parse_lipid(t, k) for (k, t) in enumerate(lipid_tables)]
    proteins = [_parse_protein(t, k) for (k, t) in enumerate(_table_array(raw, "proteins"))]
    protein_names = [p.name for p in proteins]

    surface = _table(raw, "surface")
    _check_keys(surface, _SURFACE_KEYS, "[surface]")
    grid_spacing = _get(surface, "grid_spacing", Float64, 15.0, "[surface]")
    grid_spacing > 0 || throw(ParameterError("grid_spacing must be positive"))
    grid_size = nothing
    if haskey(surface, "grid_size")
        gs = _get(surface, "grid_size", Vector{Int}, nothing, "[surface]")
        length(gs) == 2 || throw(ParameterError("grid_size must be a list of two integers [nx, ny]"))
        grid_size = (gs[1], gs[2])
    end
    mask_radius = _get(surface, "mask_radius", Float64, 0.0, "[surface]")
    mask_radius >= 0 || throw(ParameterError("mask_radius must not be negative"))
    mask_proteins = _get(surface, "mask_proteins", Vector{String},
                         mask_radius > 0 ? protein_names : String[], "[surface]")
    for p in mask_proteins
        p in protein_names || throw(ParameterError("[surface] mask_proteins: unknown protein \"$p\""))
    end
    if mask_radius > 0 && isempty(mask_proteins)
        throw(ParameterError("[surface] mask_radius > 0 but there is no protein to mask; define [[proteins]]"))
    end
    mask_residues = haskey(surface, "mask_residues") ?
                    _get(surface, "mask_residues", Vector{Int}, nothing, "[surface]") : nothing
    fill_text = _get(surface, "fill", String, mask_radius > 0 ? "min_curvature" : "neighbors", "[surface]")
    fill_text in ("neighbors", "min_curvature") ||
        throw(ParameterError("[surface] fill must be \"neighbors\" or \"min_curvature\" (got \"$fill_text\")"))
    if mask_radius > 0 && fill_text != "min_curvature"
        throw(ParameterError("[surface] masking around proteins (mask_radius > 0) requires fill = \"min_curvature\""))
    end
    center_threshold = _get(surface, "center_threshold", Float64, 20.0, "[surface]")
    (isfinite(center_threshold) && center_threshold > 0) ||
        throw(ParameterError("[surface] center_threshold must be a positive number"))
    (isfinite(grid_spacing) && isfinite(mask_radius)) ||
        throw(ParameterError("[surface] grid_spacing and mask_radius must be finite numbers"))

    spectrum = _table(raw, "spectrum")
    _check_keys(spectrum, _SPECTRUM_KEYS, "[spectrum]")
    spectrum_defaults = Dict{String,Any}(
        "q_min" => _get(spectrum, "q_min", Float64, 0.01, "[spectrum]"),
        "q_max" => _get(spectrum, "q_max", Float64, 0.21, "[spectrum]"),
        "n_bins" => _get(spectrum, "n_bins", Int, 10, "[spectrum]"))

    statistics = _table(raw, "statistics")
    _check_keys(statistics, _STATISTICS_KEYS, "[statistics]")
    blocking = _get(statistics, "blocking", String, "equal", "[statistics]")
    blocking in ("equal", "trajectory", "auto") ||
        throw(ParameterError("[statistics] blocking must be \"equal\", \"trajectory\" or \"auto\" (got \"$blocking\")"))
    n_blocks = _get(statistics, "n_blocks", Int, 10, "[statistics]")
    n_blocks >= 2 || throw(ParameterError("[statistics] n_blocks must be at least 2"))

    output = _table(raw, "output")
    _check_keys(output, _OUTPUT_KEYS, "[output]")
    directory = _get(output, "directory", String, "membrane_analysis_results", "[output]")
    save_spectra = _get(output, "save_spectra", Bool, true, "[output]")

    analyses = Dict{String,Any}[]
    for (k, t) in enumerate(_table_array(raw, "analysis"))
        loc = "[[analysis]] entry $k"
        type = _require(t, "type", String, loc)
        haskey(_ANALYSIS_KEYS, type) || throw(ParameterError(
            "unknown analysis type \"$type\" in $loc. Available types: " * join(analysis_types(), ", ")))
        _check_keys(t, vcat(["type", "label"], _ANALYSIS_KEYS[type]), "$loc (type \"$type\")")
        _validate_analysis(t, type, "$loc (type \"$type\")")
        push!(analyses, Dict{String,Any}(String(key) => val for (key, val) in t))
    end
    isempty(analyses) && throw(ParameterError("no analyses requested; add at least one [[analysis]] table"))

    return AnalysisConfig(title, lipids, proteins, grid_spacing, grid_size, Symbol(fill_text),
                          mask_proteins, mask_residues, mask_radius, center_threshold,
                          spectrum_defaults, Symbol(blocking), n_blocks, directory, save_spectra,
                          analyses, String(source))
end

# checks of the values in an [[analysis]] table that do not need the structure
function _validate_analysis(t::AbstractDict, type::AbstractString, loc::AbstractString)
    haskey(t, "direction") && _parse_direction(t["direction"], loc)
    _label_kw(t, loc)
    for k in ("q_min", "q_max", "z_cutoff")
        if haskey(t, k)
            v = _get(t, k, Float64, nothing, loc)
            isfinite(v) || throw(ParameterError("$loc: $k must be a finite number"))
        end
    end
    for k in ("n_bins", "n_fit", "residue")
        haskey(t, k) && _get(t, k, Int, nothing, loc)
    end
    haskey(t, "lipids") && _get(t, "lipids", Vector{String}, nothing, loc)
    haskey(t, "residues") && _get(t, "residues", Vector{Int}, nothing, loc)
    haskey(t, "protein") && _get(t, "protein", String, nothing, loc)
    if type in ("undulation", "lipid_sampled_curvature", "protein_sampled_curvature", "tcb")
        _require(t, "q_max", Float64, loc)
    end
    type == "tcb" && _require(t, "z_cutoff", Float64, loc)
    if type in ("protein_curvature_spectrum", "protein_thickness_spectrum")
        _require(t, "residue", Int, loc)
    end
    if haskey(t, "model")
        m = _get(t, "model", String, nothing, loc)
        m in ("HC", "MNK", "TD") ||
            throw(ParameterError("$loc: model must be \"HC\", \"MNK\" or \"TD\" (got \"$m\")"))
    end
    if haskey(t, "ref_atoms")
        t["ref_atoms"] isa AbstractDict ||
            throw(ParameterError("$loc: ref_atoms must be a table, e.g. ref_atoms = { POPC = \"C1A\" }"))
        for (key, val) in t["ref_atoms"]
            val isa AbstractString || throw(ParameterError("$loc: ref_atoms.$key must be an atom name"))
        end
    end
    if haskey(t, "weights")
        t["weights"] isa AbstractDict ||
            throw(ParameterError("$loc: weights must be a table, e.g. weights = { DOPC = 0.9, CHL1 = 0.1 }"))
        for (key, val) in t["weights"]
            val isa Real || throw(ParameterError("$loc: weights.$key must be a number"))
        end
    end
    return nothing
end

"""
    load_parameters(path) -> AnalysisConfig

Read and validate a TOML parameters file.
"""
function load_parameters(path::AbstractString)
    isfile(path) || throw(ParameterError("parameters file not found: $path"))
    raw = try
        TOML.parsefile(path)
    catch err
        err isa TOML.ParserError || rethrow()
        throw(ParameterError("could not parse $path: " * sprint(showerror, err)))
    end
    return parse_parameters(raw; source=abspath(path))
end

# --------------------------------------------------------------------------- resolution

"""
    build_analyses(config, system; box_estimate, grid_size) -> Vector{AbstractAnalysis}

Create the analysis objects requested in `config` for `system`. `box_estimate` is an estimate
`(Lx, Ly, Lz)` of the mean box and `grid_size` the surface grid; both are used to choose the
number of undulation shells.
"""
function build_analyses(cfg::AnalysisConfig, system::MembraneSystem;
                        box_estimate::NTuple{3,Float64}, grid_size::Tuple{Int,Int})
    analyses = AbstractAnalysis[]
    for (k, t) in enumerate(cfg.analyses)
        loc = "[[analysis]] entry $k (type \"$(t["type"])\")"
        a = try
            _build_analysis(t, cfg, system, box_estimate, grid_size, loc)
        catch err
            err isa ArgumentError || rethrow()
            throw(ParameterError("$loc: " * err.msg))
        end
        push!(analyses, a)
    end
    ids = analysis_id.(analyses)
    for id in unique(ids)
        count(==(id), ids) > 1 && throw(ParameterError(
            "two analyses have the same type and label ($id); give them different `label`s"))
    end
    files = lowercase.(reduce(vcat, [output_files(a, system) for a in analyses]))  # case-insensitive file systems
    for f in unique(files)
        count(==(f), files) > 1 && throw(ParameterError(
            "two analyses would write the same output file ($f); give them different `label`s"))
    end
    return analyses
end

function _spectrum_args(t, cfg, loc)
    q_min = _get(t, "q_min", Float64, cfg.spectrum_defaults["q_min"], loc)
    q_max = _get(t, "q_max", Float64, cfg.spectrum_defaults["q_max"], loc)
    n_bins = _get(t, "n_bins", Int, cfg.spectrum_defaults["n_bins"], loc)
    return (q_min=q_min, q_max=q_max, n_bins=n_bins)
end

function _protein_arg(t, system, loc)
    if haskey(t, "protein")
        return _get(t, "protein", String, nothing, loc)
    end
    length(system.proteins) == 1 && return system.proteins[1].name
    throw(ParameterError("$loc: specify `protein` (defined proteins: " *
                         (isempty(system.proteins) ? "none" : join((p.name for p in system.proteins), ", ")) * ")"))
end

function _label_kw(t, loc)
    haskey(t, "label") || return NamedTuple()
    label = _get(t, "label", String, nothing, loc)
    (isempty(label) || occursin(r"[/\\\s]", label)) &&
        throw(ParameterError("$loc: label must be non-empty and contain no spaces or slashes"))
    return (label=label,)
end

function _build_analysis(t, cfg, system, box_estimate, grid_size, loc)
    type = t["type"]
    lipids = haskey(t, "lipids") ? _get(t, "lipids", Vector{String}, nothing, loc) : nothing
    direction = _parse_direction(get(t, "direction", nothing), loc)
    label = _label_kw(t, loc)
    if type == "box"
        return BoxAnalysis(; label=get(label, :label, ""))
    elseif type == "undulation"
        q_max = _require(t, "q_max", Float64, loc)
        Lx = box_estimate[1]
        n_requested = max(floor(Int, q_max / (2π / Lx)), 4)
        n_max = min(n_requested, (minimum(grid_size) - 1) ÷ 2)   # largest non-aliased index
        n_max >= 1 || throw(ParameterError("$loc: the $(grid_size) grid is too coarse for the undulation analysis"))
        if n_max < n_requested
            @warn "$loc: undulation modes limited to n, m ≤ $n_max by the $(grid_size[1])×$(grid_size[2]) grid; " *
                  "use a smaller grid_spacing for more modes"
        end
        return UndulationAnalysis(; n_max=n_max, model=_get(t, "model", String, "TD", loc),
                                  n_fit=_get(t, "n_fit", Int, 0, loc),
                                  label=get(label, :label, ""))
    elseif type == "lipid_heights"
        return LipidHeights(system; lipids=lipids, label...)
    elseif type == "lipid_sampled_curvature"
        return LipidSampledCurvature(system; q_max=_require(t, "q_max", Float64, loc),
                                     q_min=_get(t, "q_min", Float64, 0.0, loc),
                                     direction=direction, lipids=lipids, label...)
    elseif type == "lipid_curvature_spectrum" || type == "lipid_thickness_spectrum"
        kind = type == "lipid_curvature_spectrum" ? :curvature : :thickness
        ref_atoms = Dict{String,String}()
        if haskey(t, "ref_atoms")
            for (key, val) in t["ref_atoms"]
                ref_atoms[String(key)] = String(val)
            end
        end
        return LipidSpectrum(system, kind; _spectrum_args(t, cfg, loc)..., direction=direction,
                             lipids=lipids, ref_atoms=ref_atoms, label...)
    elseif type == "protein_heights"
        residues = haskey(t, "residues") ? _get(t, "residues", Vector{Int}, nothing, loc) : nothing
        return ProteinHeights(system; protein=_protein_arg(t, system, loc), residues=residues, label...)
    elseif type == "protein_sampled_curvature"
        residues = haskey(t, "residues") ? _get(t, "residues", Vector{Int}, nothing, loc) : nothing
        return ProteinSampledCurvature(system; protein=_protein_arg(t, system, loc),
                                       q_max=_require(t, "q_max", Float64, loc),
                                       q_min=_get(t, "q_min", Float64, 0.0, loc),
                                       direction=direction, residues=residues, label...)
    elseif type == "protein_curvature_spectrum" || type == "protein_thickness_spectrum"
        kind = type == "protein_curvature_spectrum" ? :curvature : :thickness
        return ProteinSpectrum(system, kind; protein=_protein_arg(t, system, loc),
                               residue=_require(t, "residue", Int, loc),
                               _spectrum_args(t, cfg, loc)..., direction=direction, label...)
    elseif type == "tcb"
        weights = nothing
        if haskey(t, "weights")
            weights = Dict{String,Float64}(String(key) => Float64(val) for (key, val) in t["weights"])
        end
        return TCBAnalysis(system; q_max=_require(t, "q_max", Float64, loc),
                           q_min=_get(t, "q_min", Float64, 0.0, loc), direction=direction,
                           z_cutoff=_require(t, "z_cutoff", Float64, loc), weights=weights,
                           lipids=lipids, label=get(label, :label, ""))
    end
    throw(ParameterError("$loc: unknown analysis type"))
end

"""
    surface_settings(config, system, grid_size) -> SurfaceSettings

Surface settings for `config`, with mask atoms resolved in `system`.
"""
function surface_settings(cfg::AnalysisConfig, system::MembraneSystem, grid_size::Tuple{Int,Int})
    mask_atoms = Int[]
    if cfg.mask_radius > 0
        for name in cfg.mask_proteins
            k = protein_index(system, name)
            bb = system.protein_backbone[k]
            if cfg.mask_residues !== nothing
                bb = [i for i in bb if system.atoms.resnum[i] in cfg.mask_residues]
                isempty(bb) && throw(ParameterError(
                    "[surface] mask_residues: no backbone atoms of protein $name in residues $(cfg.mask_residues)"))
            end
            append!(mask_atoms, bb)
        end
    end
    return SurfaceSettings(grid_size, cfg.fill, mask_atoms, cfg.mask_radius, cfg.center_threshold)
end

"""
    config_dict(config) -> Dict{String,Any}

The settings of `config` (including defaults) as a dictionary, as written to
`parameters_used.toml`.
"""
function config_dict(cfg::AnalysisConfig)
    lipid_dict(l) = Dict{String,Any}("name" => l.name, "head_atom" => l.head_atom,
                                     "ref_atom" => l.ref_atom, "tail_atom" => l.tail_atom,
                                     "n_atoms" => l.n_atoms, "height_atoms" => l.height_atoms)
    protein_dict(p) = Dict{String,Any}("name" => p.name, "segments" => p.segments,
                                       "chains" => p.chains, "segment_prefix" => p.segment_prefix,
                                       "backbone_atoms" => p.backbone_atoms,
                                       "leaflet" => String(p.leaflet))
    surface = Dict{String,Any}("grid_spacing" => cfg.grid_spacing, "fill" => String(cfg.fill),
                               "mask_proteins" => cfg.mask_proteins, "mask_radius" => cfg.mask_radius,
                               "center_threshold" => cfg.center_threshold)
    cfg.grid_size === nothing || (surface["grid_size"] = collect(cfg.grid_size))
    cfg.mask_residues === nothing || (surface["mask_residues"] = cfg.mask_residues)
    d = Dict{String,Any}(
        "lipids" => [lipid_dict(l) for l in cfg.lipids],
        "surface" => surface,
        "spectrum" => cfg.spectrum_defaults,
        "statistics" => Dict{String,Any}("blocking" => String(cfg.blocking), "n_blocks" => cfg.n_blocks),
        "output" => Dict{String,Any}("directory" => cfg.output_directory, "save_spectra" => cfg.save_spectra),
        "analysis" => cfg.analyses)
    isempty(cfg.title) || (d["title"] = cfg.title)
    isempty(cfg.proteins) || (d["proteins"] = [protein_dict(p) for p in cfg.proteins])
    return d
end
