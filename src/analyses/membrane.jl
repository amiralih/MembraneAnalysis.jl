# Box dimensions / area compressibility, and the undulation (height fluctuation) spectrum.

_label_suffix(label::AbstractString) = isempty(label) ? "" : "_" * label

# ------------------------------------------------------------------------------------- box

"""
    BoxAnalysis(; label="")

Mean box lengths, membrane area `A = Lx Ly` and the area expansion (compressibility) modulus
`K_A = ⟨A⟩ / (⟨A²⟩ − ⟨A⟩²)` in k_BT/Å².

Output files: `box_dims.dat`, `A.dat`, `KA.dat` (with `_<label>` appended when a label is set).
"""
struct BoxAnalysis <: AbstractAnalysis
    label::String
end
BoxAnalysis(; label::AbstractString="") = BoxAnalysis(String(label))

analysis_id(a::BoxAnalysis) = "box" * _label_suffix(a.label)
output_files(a::BoxAnalysis, system) =
    ["box_dims$(_label_suffix(a.label)).dat", "A$(_label_suffix(a.label)).dat",
     "KA$(_label_suffix(a.label)).dat"]
component_labels(::BoxAnalysis) = ["Lx", "Ly", "Lz", "A", "A^2"]
needs_surface(::BoxAnalysis) = false

function accumulate_frame!(sums, counts, ::BoxAnalysis, system, ctx::FrameContext)
    A = ctx.Lx * ctx.Ly
    sums[1] += ctx.Lx
    sums[2] += ctx.Ly
    sums[3] += ctx.Lz
    sums[4] += A
    sums[5] += A^2
    counts .+= 1
    return nothing
end

_area_modulus(v) = v[4] / (v[5] - v[4]^2)

function write_results(a::BoxAnalysis, bs::BlockStats, system, out::OutputContext)
    m = block_means(bs)
    e = block_sems(bs)
    sfx = _label_suffix(a.label)
    write_table(out, "box_dims$(sfx).dat", "mean box lengths [Å]; " * stats_header(bs),
                ["quantity", "mean", "sem"], ["Lx", "Ly", "Lz"], m[1:3], e[1:3])
    write_table(out, "A$(sfx).dat", "membrane area A = Lx Ly [Å^2]; " * stats_header(bs),
                ["mean", "sem"], [m[4]], [e[4]])
    KA, KA_err = jackknife(_area_modulus, bs)
    write_table(out, "KA$(sfx).dat",
                "area expansion modulus K_A = <A>/(<A^2> - <A>^2) [k_BT/Å^2]; " * stats_header(bs),
                ["KA", "error"], [KA], [KA_err])
    return Dict{String,Any}("Lx" => m[1], "Ly" => m[2], "Lz" => m[3], "Lx_sem" => e[1],
                            "Ly_sem" => e[2], "Lz_sem" => e[3], "A" => m[4], "A_sem" => e[4],
                            "KA" => KA, "KA_error" => KA_err)
end

# ------------------------------------------------------------------------------ undulation

"""
    UndulationAnalysis(; n_max, model="TD", n_fit=0, label="")

Mean squared height and thickness modes, ⟨|h_q|²⟩ and ⟨|t_q|²⟩, in shells of equal
`|(n, m)|` (`0 ≤ n, m ≤ n_max`, see [`get_index_pairs`](@ref)), and the bending modulus from a
fit of `q⁴⟨|h_q|²⟩` to one of the models of Terzi & Deserno, J. Chem. Phys. 147, 084702 (2017):
`"HC"` (Helfrich–Canham, eq. 36), `"MNK"` (May–Narang–Kopelevich, eq. 37) or `"TD"`
(Terzi–Deserno, eq. 47). `n_fit` limits the fit to the first `n_fit` shells (0: all).
`k_c = A / p₁` in k_BT, with a block-jackknife error.

Output files: `hq2.dat`, `tq2.dat`, `kc_hq2.dat` (with `_<label>` appended when set).
"""
struct UndulationAnalysis <: AbstractAnalysis
    label::String
    n_max::Int
    shell_values::Vector{Float64}
    shell_pairs::Vector{Vector{Tuple{Int,Int}}}
    model::String
    n_fit::Int
    labels::Vector{String}
end

function UndulationAnalysis(; n_max::Integer, model::AbstractString="TD", n_fit::Integer=0,
                            label::AbstractString="")
    model in ("HC", "MNK", "TD") || throw(ArgumentError("undulation model must be \"HC\", \"MNK\" or \"TD\" (got \"$model\")"))
    n_max >= 1 || throw(ArgumentError("n_max must be ≥ 1"))
    values, pairs = get_index_pairs(n_max)
    values = values[2:end]          # drop the q = 0 mode
    pairs = pairs[2:end]
    ns = length(values)
    labels = vcat(["hq2_$(k)" for k in 1:ns], ["tq2_$(k)" for k in 1:ns], ["Lx", "A"])
    return UndulationAnalysis(String(label), Int(n_max), values, pairs, String(model),
                              Int(n_fit), labels)
end

analysis_id(a::UndulationAnalysis) = "undulation" * _label_suffix(a.label)
output_files(a::UndulationAnalysis, system) =
    ["hq2$(_label_suffix(a.label)).dat", "tq2$(_label_suffix(a.label)).dat",
     "kc_hq2$(_label_suffix(a.label)).dat"]
component_labels(a::UndulationAnalysis) = a.labels

function accumulate_frame!(sums, counts, a::UndulationAnalysis, system, ctx::FrameContext)
    ns = length(a.shell_values)
    nx, ny = size(ctx.hq)
    for s in 1:ns
        h = 0.0
        t = 0.0
        for (i, j) in a.shell_pairs[s]
            (i <= (nx - 1) ÷ 2 && j <= (ny - 1) ÷ 2) || throw(ArgumentError(
                "undulation mode ($i, $j) is not resolved by the $(nx)×$(ny) grid; reduce q_max or grid_spacing"))
            h += abs2(ctx.hq[i + 1, j + 1])
            t += abs2(ctx.tq[i + 1, j + 1])
        end
        np = length(a.shell_pairs[s])
        sums[s] += h / np
        sums[ns + s] += t / np
    end
    sums[2ns + 1] += ctx.Lx
    sums[2ns + 2] += ctx.Lx * ctx.Ly
    counts .+= 1
    return nothing
end

# models of q⁴⟨|h_q|²⟩ (same expressions and starting values as earlier versions)
_hq2_model(::Val{:HC}) = (x, p) -> p[1] .+ 0 .* x
_hq2_model(::Val{:MNK}) = (x, p) -> @. p[1] + p[2] * x^2
_hq2_model(::Val{:TD}) = (x, p) -> @. (p[1] + p[2] * x^2) * (1 / (1 - x^2 / p[3]))
_hq2_start(model::String) = model == "HC" ? [1e3] : model == "MNK" ? [1e3, 1e5] : [1e3, 1e5, 1e-1]

function _undulation_kc(a::UndulationAnalysis, v)
    ns = length(a.shell_values)
    Lx = v[2ns + 1]
    A = v[2ns + 2]
    n = a.n_fit > 0 ? min(a.n_fit, ns) : ns
    q = a.shell_values[1:n] .* (2π / Lx)
    y = v[1:n] .* q .^ 4
    all(isfinite, y) || return NaN
    try
        fit = LsqFit.curve_fit(_hq2_model(Val(Symbol(a.model))), q, y, _hq2_start(a.model))
        return A / LsqFit.coef(fit)[1]
    catch err
        err isa InterruptException && rethrow()
        return NaN
    end
end

function write_results(a::UndulationAnalysis, bs::BlockStats, system, out::OutputContext)
    m = block_means(bs)
    e = block_sems(bs)
    ns = length(a.shell_values)
    q = a.shell_values .* (2π / m[2ns + 1])
    sfx = _label_suffix(a.label)
    write_table(out, "hq2$(sfx).dat",
                "mean squared height modes <|h_q|^2> [Å^6] in shells of equal |q|; " * stats_header(bs),
                ["q [1/Å]", "mean", "sem"], q, m[1:ns], e[1:ns])
    write_table(out, "tq2$(sfx).dat",
                "mean squared thickness modes <|t_q|^2> [Å^6] in shells of equal |q|; " * stats_header(bs),
                ["q [1/Å]", "mean", "sem"], q, m[(ns + 1):(2ns)], e[(ns + 1):(2ns)])
    kc, kc_err = jackknife(v -> _undulation_kc(a, v), bs)
    n = a.n_fit > 0 ? min(a.n_fit, ns) : ns
    write_table(out, "kc_hq2$(sfx).dat",
                "bending modulus [k_BT] from a fit of q^4 <|h_q|^2> (model $(a.model), first $n shells); " *
                stats_header(bs), ["kc", "error"], [kc], [kc_err])
    return Dict{String,Any}("kc" => kc, "kc_error" => kc_err, "model" => a.model,
                            "n_fit" => n)
end
