# Transverse curvature bias (TCB) analysis.

"""
    TCBAnalysis(system; q_max, q_min=0.0, direction=:all, z_cutoff, weights=nothing,
                lipids=nothing, label="")

Transverse curvature bias analysis (Sapp, Beaven & Sodt, Phys. Rev. E 103, 042413 (2021)).
For every lipid species the mean sampled curvature ⟨c⟩ of each height atom (see
[`LipidSampledCurvature`](@ref)) is fitted linearly against the atom's mean height z above the
midplane, using atoms with `z < z_cutoff`. With area-fraction `weights` (equal by default) the
weighted slope `s` and intercept `i` give

- the bending modulus `k_c = -2 / (A s)` (k_BT),
- the neutral surface height `-i / s` (Å),
- the curvature bias of each lipid relative to the weighted average, `Δc = i_lipid - i`.

Errors are block-jackknife estimates.

Output files (with `_X`/`_Y` for directional analyses and `_<label>` when a label is set):
`TCB.dat` (kc, neutral surface), `TCB_lipids.dat` (per-lipid fit parameters and Δc) and
`<LIPID>_TCB_profile.dat` (z and ⟨c⟩ of every atom, for plotting).
"""
struct TCBAnalysis <: AbstractAnalysis
    label::String
    tag::String
    species::Vector{Int}
    groups::Vector{Vector{AtomGroup}}
    weights::Vector{Float64}
    shell::QBins
    direction::Symbol
    z_cutoff::Float64
    labels::Vector{String}
end

function TCBAnalysis(system::MembraneSystem; q_max::Real, q_min::Real=0.0,
                     direction::Symbol=:all, z_cutoff::Real, weights=nothing, lipids=nothing,
                     label::AbstractString="")
    _check_direction(direction)
    species = _lipid_species(system, lipids)
    names = [system.lipids[s].name for s in species]
    if weights === nothing
        w = fill(1.0 / length(species), length(species))
    else
        weights = Dict{String,Float64}(String(k) => Float64(v) for (k, v) in weights)
        for k in keys(weights)
            String(k) in names || throw(ArgumentError("TCB weights: \"$k\" is not one of the analysed lipids $(names)"))
        end
        for n in names
            haskey(weights, n) || throw(ArgumentError("TCB weights: no weight given for lipid $n " *
                                                      "(give weights for all analysed lipids)"))
        end
        w = [Float64(weights[n]) for n in names]
        (all(>=(0), w) && sum(w) > 0) || throw(ArgumentError("TCB weights must be non-negative and not all zero"))
    end
    groups = [height_atom_groups(system, s) for s in species]
    labels = String[]
    for (n, gs) in zip(names, groups)
        append!(labels, ["$n:z:$(g.label)" for g in gs])
        append!(labels, ["$n:c:$(g.label)" for g in gs])
    end
    push!(labels, "A")
    return TCBAnalysis(String(label), _species_tag(system, species), species, groups, w,
                       q_shell(q_min, q_max), direction, Float64(z_cutoff), labels)
end

_tcb_suffix(a::TCBAnalysis) = direction_suffix(a.direction) * _label_suffix(a.label)
analysis_id(a::TCBAnalysis) = "tcb" * _tcb_suffix(a) * a.tag
output_files(a::TCBAnalysis, system) =
    vcat(["TCB$(_tcb_suffix(a)).dat", "TCB$(_tcb_suffix(a))_lipids.dat"],
         ["$(system.lipids[s].name)_TCB$(_tcb_suffix(a))_profile.dat" for s in a.species])
component_labels(a::TCBAnalysis) = a.labels

function accumulate_frame!(sums, counts, a::TCBAnalysis, system, ctx::FrameContext)
    field = FourierField(ctx.hq, ctx.Lx, ctx.Ly, a.shell, a.direction, :curvature)
    c = 0
    for gs in a.groups
        for g in gs
            c += 1
            _accumulate_heights!(sums, counts, c, g, system, ctx)
        end
        for g in gs
            _accumulate_field!(sums, counts, c, field, g, system, ctx, true)
            c += 1
        end
    end
    sums[c + 1] += ctx.Lx * ctx.Ly
    counts[c + 1] += 1
    return nothing
end

# component ranges (z, c) of species k
function _tcb_ranges(a::TCBAnalysis, k::Int)
    c = 0
    for j in 1:(k - 1)
        c += 2 * length(a.groups[j])
    end
    g = length(a.groups[k])
    return (c + 1):(c + g), (c + g + 1):(c + 2g)
end

"""
Fit of ⟨c⟩(z) for every species; returns
`[kc, neutral_surface, intercepts..., slopes..., Δc...]` (NaN if a fit is impossible).
"""
function _tcb_estimate(a::TCBAnalysis, v::AbstractVector)
    ns = length(a.species)
    intercepts = zeros(ns)
    slopes = zeros(ns)
    nan_result = fill(NaN, 2 + 3ns)
    for k in 1:ns
        rz, rc = _tcb_ranges(a, k)
        z = v[rz]
        c = v[rc]
        sel = [isfinite(z[i]) && isfinite(c[i]) && z[i] < a.z_cutoff for i in eachindex(z)]
        count(sel) >= 2 || return nan_result
        X = hcat(ones(count(sel)), z[sel])
        p = X \ c[sel]
        intercepts[k] = p[1]
        slopes[k] = p[2]
    end
    w = a.weights ./ sum(a.weights)
    slope = sum(w .* slopes)
    intercept = sum(w .* intercepts)
    A = v[end]
    kc = -2 / (A * slope)
    neutral_surface = -intercept / slope
    return vcat(kc, neutral_surface, intercepts, slopes, intercepts .- intercept)
end

function write_results(a::TCBAnalysis, bs::BlockStats, system, out::OutputContext)
    m = block_means(bs)
    e = block_sems(bs)
    est, err = jackknife(v -> _tcb_estimate(a, v), bs)
    ns = length(a.species)
    names = [system.lipids[s].name for s in a.species]
    sfx = _tcb_suffix(a)
    shell_text = "$(a.shell.lo[1]) < q < $(a.shell.hi[1]) 1/Å, $(_direction_text(a.direction)), " *
                 "z_cutoff = $(a.z_cutoff) Å"
    write_table(out, "TCB$(sfx).dat",
                "transverse curvature bias analysis ($shell_text); kc [k_BT], neutral surface [Å]; " *
                stats_header(bs), ["quantity", "value", "error"], ["kc", "neutral_surface"],
                est[1:2], err[1:2])
    write_table(out, "TCB$(sfx)_lipids.dat",
                "TCB fit <c> = intercept + slope z per lipid and curvature bias delta_c = " *
                "intercept - weighted mean intercept [1/Å]; " * stats_header(bs),
                ["lipid", "weight", "intercept", "intercept_err", "slope", "slope_err",
                 "delta_c", "delta_c_err"],
                names, a.weights ./ sum(a.weights), est[3:(2 + ns)], err[3:(2 + ns)],
                est[(3 + ns):(2 + 2ns)], err[(3 + ns):(2 + 2ns)], est[(3 + 2ns):(2 + 3ns)],
                err[(3 + 2ns):(2 + 3ns)])
    for (k, n) in enumerate(names)
        rz, rc = _tcb_ranges(a, k)
        write_table(out, "$(n)_TCB$(sfx)_profile.dat",
                    "$n: mean height z [Å] and sampled curvature <c> [1/Å] of atoms ($shell_text); " *
                    stats_header(bs), ["atom", "z", "z_sem", "c", "c_sem"],
                    [g.label for g in a.groups[k]], m[rz], e[rz], m[rc], e[rc])
    end
    return Dict{String,Any}("kc" => est[1], "kc_error" => err[1],
                            "neutral_surface" => est[2], "neutral_surface_error" => err[2])
end
