# The analysis interface.
#
# An analysis is a struct (subtype of `AbstractAnalysis`) holding its settings and the atom
# groups it needs. During the single pass over a trajectory, `accumulate_frame!` is called once
# per frame with column views `sums` and `counts` (one entry per component); afterwards the
# per-frame sums/counts of all trajectories are reduced to means ± errors by `write_results`.
#
# To add a new analysis, define a struct `<: AbstractAnalysis` and methods for
# `analysis_id`, `component_labels`, `output_files`, `accumulate_frame!` and `write_results`
# (and `needs_surface` if it does not use the spectra or the lipid leaflets). To make it
# available in parameters files, add its keys to `_ANALYSIS_KEYS` and a branch to
# `_build_analysis` in parameters.jl.

"""
    AbstractAnalysis

Supertype of all analyses. See the developer notes in `src/analyses/interface.jl`.
"""
abstract type AbstractAnalysis end

"""
    FrameContext

Everything an analysis can use for one frame: wrapped and centred coordinates (`coords`,
3 × N), box lengths `Lx`, `Ly`, `Lz`, the bilayer midplane height `midplane`, the leaflet
(±1) of every lipid molecule `leaflet`, and the height and thickness spectra `hq`, `tq`.
"""
struct FrameContext
    coords::Matrix{Float64}
    Lx::Float64
    Ly::Float64
    Lz::Float64
    midplane::Float64
    leaflet::Vector{Int8}
    hq::Matrix{ComplexF64}
    tq::Matrix{ComplexF64}
end

"""
    analysis_id(a) -> String

Unique identifier of an analysis instance (used for HDF5 group names and log messages).
"""
function analysis_id end

"""
    component_labels(a) -> Vector{String}

Labels of the per-frame components accumulated by an analysis.
"""
function component_labels end

"""
    output_files(a, system) -> Vector{String}

Names of the result files written by an analysis (used to detect collisions).
"""
function output_files end

"""
    accumulate_frame!(sums, counts, a, system, ctx::FrameContext)

Add the contributions of one frame to `sums` and `counts` (vectors with one entry per
component, initialised to zero for every frame).
"""
function accumulate_frame! end

"""
    write_results(a, stats::BlockStats, system, out::OutputContext) -> Dict{String,Any}

Write the result files of an analysis and return scalar results for `summary.toml`.
"""
function write_results end

"""
    needs_surface(a) -> Bool

Whether an analysis needs the fluctuation spectra and lipid leaflets (default `true`).
"""
needs_surface(::AbstractAnalysis) = true

n_components(a::AbstractAnalysis) = length(component_labels(a))

"""
    direction_suffix(direction) -> String

File-name suffix of a direction (`""`, `"_X"` or `"_Y"`).
"""
direction_suffix(d::Symbol) = d === :x ? "_X" : d === :y ? "_Y" : ""

# leaflet sign of atom `idx`: dynamic leaflet for lipids, fixed sign for protein groups
@inline _sign(ctx::FrameContext, system::MembraneSystem, g::AtomGroup, k::Int, idx::Int) =
    isempty(g.signs) ? ctx.leaflet[system.atom_mol[idx]] : g.signs[k]

# add the (signed) height from the midplane of the atoms of group `g` to component `c`
function _accumulate_heights!(sums, counts, c::Int, g::AtomGroup, system, ctx::FrameContext)
    s = 0.0
    @inbounds for (k, idx) in enumerate(g.atoms)
        s += _sign(ctx, system, g, k, idx) * (ctx.coords[3, idx] - ctx.midplane)
    end
    sums[c] += s
    counts[c] += length(g.atoms)
    return nothing
end

# add the per-mode values of `field` at the atoms of group `g` to components c0+1 … c0+nbins
function _accumulate_field!(sums, counts, c0::Int, field::FourierField, g::AtomGroup, system,
                            ctx::FrameContext, signed::Bool)
    nb = length(field.n_modes)
    @inbounds for (k, idx) in enumerate(g.atoms)
        v = evaluate!(field, ctx.coords[1, idx], ctx.coords[2, idx])
        sg = signed ? _sign(ctx, system, g, k, idx) : Int8(1)
        for b in 1:nb
            field.n_modes[b] == 0 && continue
            sums[c0 + b] += sg * v[b]
            counts[c0 + b] += 1
        end
    end
    return nothing
end
