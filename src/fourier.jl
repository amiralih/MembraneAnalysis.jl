# Evaluation of curvature and thickness at points from the fluctuation spectra, restricted to
# shells (bins) of wave-vector magnitude.

"""
    QBins(lo, hi, centers)

Shells of wave-vector magnitude `lo[k] < |q| < hi[k]` (strict inequalities, in Å⁻¹).
Use [`spectrum_bins`](@ref) or [`q_shell`](@ref) to construct them.
"""
struct QBins
    lo::Vector{Float64}
    hi::Vector{Float64}
    centers::Vector{Float64}
end

Base.length(b::QBins) = length(b.lo)

"""
    spectrum_bins(q_min, q_max, n_bins) -> QBins

`n_bins` equal-width shells between `q_min` and `q_max`.
"""
function spectrum_bins(q_min::Real, q_max::Real, n_bins::Integer)
    (q_max > q_min >= 0 && n_bins >= 1) ||
        throw(ArgumentError("need 0 ≤ q_min < q_max and n_bins ≥ 1 (got $q_min, $q_max, $n_bins)"))
    dq = (q_max - q_min) / n_bins
    lo = [q_min + (k - 1) * dq for k in 1:n_bins]
    hi = [q_min + k * dq for k in 1:n_bins]
    centers = [q_min + (k - 0.5) * dq for k in 1:n_bins]
    return QBins(lo, hi, centers)
end

"""
    q_shell(q_min, q_max) -> QBins

A single shell `q_min < |q| < q_max`.
"""
function q_shell(q_min::Real, q_max::Real)
    q_max > q_min >= 0 || throw(ArgumentError("need 0 ≤ q_min < q_max (got $q_min, $q_max)"))
    return QBins([Float64(q_min)], [Float64(q_max)], [(q_min + q_max) / 2])
end

"""
    FourierField(spectrum, Lx, Ly, bins, direction, kind)

Prepared evaluation of a spectrum at arbitrary points for one frame.

For every bin `b` the value at a point `r` is

    f_b(r) = Re[ ∑_{q ∈ b} w(q) s_q e^{i q·r} ] / (Lx Ly) / (n_b / 2)

where `s_q` is the spectrum, `n_b` the number of modes in the bin (so the result is per mode),
and `w(q) = |q|²` for `kind = :curvature` (curvature of the surface) or `1` for
`kind = :thickness`. Modes with `|i| ≤ nx÷2 - 1`, `|j| ≤ ny÷2 - 1` are used; `direction = :x`
(`:y`) keeps only modes with wave vectors along x (y).
"""
struct FourierField
    i::Vector{Int}
    j::Vector{Int}
    bin::Vector{Int}
    coef::Vector{ComplexF64}
    n_modes::Vector{Int}      # number of modes (degrees of freedom) in every bin
    kx::Float64               # 2π / Lx
    ky::Float64               # 2π / Ly
    mx::Int
    my::Int
    ex::Vector{ComplexF64}    # scratch: e^{i kx i x}
    ey::Vector{ComplexF64}    # scratch: e^{i ky j y}
    values::Vector{Float64}   # scratch: values per bin at the last point
end

function FourierField(spectrum::AbstractMatrix{<:Complex}, Lx::Real, Ly::Real, bins::QBins,
                      direction::Symbol, kind::Symbol)
    nx, ny = size(spectrum)
    mx = nx ÷ 2 - 1
    my = ny ÷ 2 - 1
    weighted = kind === :curvature ? true : kind === :thickness ? false :
               throw(ArgumentError("kind must be :curvature or :thickness"))
    direction in (:all, :x, :y) || throw(ArgumentError("direction must be :all, :x or :y"))
    nb = length(bins)
    ii = Int[]; jj = Int[]; bb = Int[]; coef = ComplexF64[]
    n_modes = zeros(Int, nb)
    for i in -mx:mx, j in -my:my
        direction === :x && j != 0 && continue
        direction === :y && i != 0 && continue
        qx = (2π / Lx) * i
        qy = (2π / Ly) * j
        q2 = qx^2 + qy^2
        for b in 1:nb
            if bins.lo[b]^2 < q2 < bins.hi[b]^2
                push!(ii, i); push!(jj, j); push!(bb, b)
                s = spectrum[mod(i, nx) + 1, mod(j, ny) + 1] / (Lx * Ly)
                push!(coef, weighted ? s * q2 : s)
                n_modes[b] += 1
                break
            end
        end
    end
    return FourierField(ii, jj, bb, coef, n_modes, 2π / Lx, 2π / Ly, mx, my,
                        zeros(ComplexF64, 2mx + 1), zeros(ComplexF64, 2my + 1), zeros(nb))
end

"""
    evaluate!(field, x, y) -> Vector{Float64}

Values of `field` in every bin at the point `(x, y)` (`NaN` for bins without modes). The
returned vector is reused by the next call.
"""
function evaluate!(f::FourierField, x::Real, y::Real)
    @inbounds for i in -f.mx:f.mx
        f.ex[i + f.mx + 1] = cis((f.kx * i) * x)
    end
    @inbounds for j in -f.my:f.my
        f.ey[j + f.my + 1] = cis((f.ky * j) * y)
    end
    fill!(f.values, 0.0)
    @inbounds for m in eachindex(f.coef)
        f.values[f.bin[m]] += real(f.coef[m] * f.ex[f.i[m] + f.mx + 1] * f.ey[f.j[m] + f.my + 1])
    end
    @inbounds for b in eachindex(f.values)
        f.values[b] /= f.n_modes[b] / 2
    end
    return f.values
end

"""
    curvature_at(hq, Lx, Ly, x, y; q_min=0, q_max, direction=:all) -> Float64

Mean curvature per mode at `(x, y)` from the height spectrum `hq`, using modes with
`q_min < |q| < q_max`.
"""
curvature_at(hq, Lx, Ly, x, y; q_min=0.0, q_max, direction::Symbol=:all) =
    evaluate!(FourierField(hq, Lx, Ly, q_shell(q_min, q_max), direction, :curvature), x, y)[1]

"""
    thickness_at(tq, Lx, Ly, x, y; q_min=0, q_max, direction=:all) -> Float64

Thickness deviation per mode at `(x, y)` from the thickness spectrum `tq`, using modes with
`q_min < |q| < q_max`.
"""
thickness_at(tq, Lx, Ly, x, y; q_min=0.0, q_max, direction::Symbol=:all) =
    evaluate!(FourierField(tq, Lx, Ly, q_shell(q_min, q_max), direction, :thickness), x, y)[1]

"""
    get_index_pairs(N) -> (values, index_pairs)

All index pairs `(n, m)` with `0 ≤ n, m ≤ N` grouped by `v = √(n² + m²) ≤ N`, sorted by `v`.
Used to group height-spectrum modes into shells of equal |q|.
"""
function get_index_pairs(N::Integer)
    all_values = Float64[]
    all_pairs = Tuple{Int,Int}[]
    for i in 0:N, j in 0:N
        push!(all_pairs, (i, j))
        push!(all_values, sqrt(i^2 + j^2))
    end
    order = sortperm(all_values; alg=MergeSort)   # stable: keeps (n, m) order within a shell
    values = Float64[]
    index_pairs = Vector{Tuple{Int,Int}}[]
    for k in order
        all_values[k] > N && break
        if !isempty(values) && all_values[k] == values[end]
            push!(index_pairs[end], all_pairs[k])
        else
            push!(values, all_values[k])
            push!(index_pairs, [all_pairs[k]])
        end
    end
    return values, index_pairs
end
