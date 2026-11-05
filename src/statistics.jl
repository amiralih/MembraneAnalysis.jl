# Block averaging, block jackknife and automatic block-size selection.
#
# Every analysis produces, for every frame, a vector of sums and a vector of counts (one entry
# per "component", e.g. one q bin). The mean of a component is ∑ sums / ∑ counts over the
# frames. Statistical errors are obtained by dividing the (time-ordered) frames into blocks.

"""
    blocking_error(a, n_blocks=10)

Standard error of the mean of the time series `a`, estimated from `n_blocks` contiguous blocks
of `length(a) ÷ n_blocks` points (the remaining points are not used for the error).
"""
function blocking_error(a::AbstractVector{<:Real}, n_blocks::Integer=10)
    n = length(a)
    m = n ÷ n_blocks
    if m < 1
        @warn "too few data points ($n) for $n_blocks blocks; returning the naive standard error"
        return std(a) / sqrt(n)
    end
    block_means = [mean(view(a, (1 + (b - 1) * m):(b * m))) for b in 1:n_blocks]
    return std(block_means) / sqrt(n_blocks)
end

"""
    equal_blocks(n_frames, n_blocks) -> Vector{UnitRange{Int}}

`n_blocks` contiguous blocks of `n_frames ÷ n_blocks` frames, starting at the first frame
(the last `n_frames % n_blocks` frames are not part of any block). With fewer frames than
blocks every frame is its own block.
"""
function equal_blocks(n_frames::Integer, n_blocks::Integer)
    n_blocks >= 2 || throw(ArgumentError("at least 2 blocks are needed (got $n_blocks)"))
    m = n_frames ÷ n_blocks
    if m < 1
        @warn "only $n_frames frames for $n_blocks blocks; using one block per frame" maxlog = 1
        return [i:i for i in 1:n_frames]
    end
    return [(1 + (b - 1) * m):(b * m) for b in 1:n_blocks]
end

"""
    blocks_of_length(n_frames, len) -> Vector{UnitRange{Int}}

Contiguous blocks of `len` frames.
"""
blocks_of_length(n_frames::Integer, len::Integer) =
    [(1 + (b - 1) * len):(b * len) for b in 1:(n_frames ÷ len)]

"""
    trajectory_blocks(n_frames_per_trajectory) -> Vector{UnitRange{Int}}

One block per trajectory (frames of all trajectories concatenated in order).
"""
function trajectory_blocks(lengths::AbstractVector{<:Integer})
    stops = cumsum(lengths)
    return [(stop - len + 1):stop for (len, stop) in zip(lengths, stops)]
end

# 99 % quantiles of the χ² distribution with 1, 2, …, 64 degrees of freedom
const CHI2_QUANTILES_99 = (
    6.634897, 9.210340, 11.344867, 13.276704, 15.086272, 16.811894, 18.475307, 20.090235,
    21.665994, 23.209251, 24.724970, 26.216967, 27.688250, 29.141238, 30.577914, 31.999927,
    33.408664, 34.805306, 36.190869, 37.566235, 38.932173, 40.289360, 41.638398, 42.979820,
    44.314105, 45.641683, 46.962942, 48.278236, 49.587884, 50.892181, 52.191395, 53.485772,
    54.775540, 56.060909, 57.342073, 58.619215, 59.892500, 61.162087, 62.428121, 63.690740,
    64.950071, 66.206236, 67.459348, 68.709513, 69.956832, 71.201400, 72.443307, 73.682639,
    74.919474, 76.153891, 77.385962, 78.615756, 79.843338, 81.068772, 82.292117, 83.513430,
    84.732766, 85.950176, 87.165711, 88.379419, 89.591344, 90.801532, 92.010024, 93.216860,
)

"""
    optimal_block_length(x) -> Union{Int, Nothing}

Block length (a power of two) after which successive blocking transformations of the time
series `x` no longer change the error estimate, determined with the automated blocking method
of M. Jonsson, Phys. Rev. E 98, 043304 (2018) (Flyvbjerg–Petersen blocking with a χ² test at
the 1 % level). Uses the first `2^⌊log₂ n⌋` points. Returns `nothing` if `x` is too short,
constant, or contains non-finite values.
"""
function optimal_block_length(x::AbstractVector{<:Real})
    n = length(x)
    n < 16 && return nothing
    all(isfinite, x) || return nothing
    d = floor(Int, log2(n))
    y = Float64.(x[1:(2^d)])
    μ = mean(y)
    s = zeros(d)
    γ = zeros(d)
    nk = zeros(d)
    for k in 1:d
        m = length(y)
        nk[k] = m
        γ[k] = sum((y[t] - μ) * (y[t + 1] - μ) for t in 1:(m - 1)) / m
        s[k] = var(y; corrected=false)
        y = (y[1:2:end] .+ y[2:2:end]) ./ 2
    end
    any(<=(0), s) && return nothing
    M = reverse(cumsum(reverse(nk .* (γ ./ s) .^ 2)))
    for k in 1:d
        df = d - k + 1
        if M[k] < CHI2_QUANTILES_99[min(df, length(CHI2_QUANTILES_99))]
            return 2^(k - 1)
        end
    end
    return 2^(d - 1)
end

"""
    BlockStats

Per-component totals and per-block sums of an analysis, from which means, standard errors and
jackknife estimates are computed. Construct with `BlockStats(labels, sums, counts, blocks)`,
where `sums` and `counts` are `n_components × n_frames` matrices.
"""
struct BlockStats
    labels::Vector{String}
    total_sum::Vector{Float64}
    total_count::Vector{Float64}
    block_sum::Matrix{Float64}
    block_count::Matrix{Float64}
    blocks::Vector{UnitRange{Int}}
end

function BlockStats(labels::Vector{String}, sums::AbstractMatrix, counts::AbstractMatrix,
                    blocks::Vector{UnitRange{Int}})
    size(sums) == size(counts) || throw(DimensionMismatch("sums and counts differ in size"))
    nc = size(sums, 1)
    nb = length(blocks)
    block_sum = zeros(nc, nb)
    block_count = zeros(nc, nb)
    for (b, r) in enumerate(blocks)
        block_sum[:, b] .= vec(sum(view(sums, :, r); dims=2))
        block_count[:, b] .= vec(sum(view(counts, :, r); dims=2))
    end
    return BlockStats(labels, vec(sum(sums; dims=2)), vec(sum(counts; dims=2)),
                      block_sum, block_count, blocks)
end

n_blocks(bs::BlockStats) = length(bs.blocks)

"""
    block_means(bs) -> Vector{Float64}

Mean of every component over all frames (∑ sums / ∑ counts; `NaN` if never sampled).
"""
block_means(bs::BlockStats) = bs.total_sum ./ bs.total_count

"""
    block_sems(bs) -> Vector{Float64}

Standard error of the mean of every component: standard deviation of the block means divided
by √(number of blocks). Blocks in which a component was not sampled are ignored.
"""
function block_sems(bs::BlockStats)
    bm = bs.block_sum ./ bs.block_count
    return [_sem(view(bm, c, :)) for c in axes(bm, 1)]
end

function _sem(v)
    w = filter(isfinite, collect(v))
    length(w) < 2 && return NaN
    return std(w) / sqrt(length(w))
end

"""
    jackknife(f, bs) -> (estimate, error)

Block-jackknife estimate and standard error of a derived quantity `f(means)`, where `means` is
the vector of component means. `f` may return a number or a vector. Leave-one-block-out samples
for which `f` returns non-finite values are ignored.
"""
function jackknife(f, bs::BlockStats)
    θ = f(block_means(bs))
    samples = typeof(θ)[]
    for b in 1:n_blocks(bs)
        m = (bs.total_sum .- view(bs.block_sum, :, b)) ./ (bs.total_count .- view(bs.block_count, :, b))
        t = f(m)
        all(isfinite, t) && push!(samples, t)
    end
    nb = length(samples)
    nb < 2 && return θ, θ .* NaN
    t̄ = sum(samples) ./ nb
    err = sqrt.((nb - 1) / nb .* sum(t -> (t .- t̄) .^ 2, samples))
    return θ, err
end

"""
    choose_blocks(mode, sums, counts, trajectory_lengths; n_blocks=10) -> Vector{UnitRange{Int}}

Blocks used for the statistics of one analysis:

- `:equal` — [`equal_blocks`](@ref)`(n_frames, n_blocks)`;
- `:trajectory` — one block per trajectory;
- `:auto` — blocks of the largest [`optimal_block_length`](@ref) over all components of the
  analysis (falls back to `:equal` when there are too few frames).
"""
function choose_blocks(mode::Symbol, sums::AbstractMatrix, counts::AbstractMatrix,
                       trajectory_lengths::AbstractVector{<:Integer}; n_blocks::Integer=10)
    n = size(sums, 2)
    if mode === :equal
        return equal_blocks(n, n_blocks)
    elseif mode === :trajectory
        if length(trajectory_lengths) < 2
            @warn "block averaging per trajectory needs at least 2 trajectories; using $n_blocks equal blocks" maxlog = 1
            return equal_blocks(n, n_blocks)
        end
        return trajectory_blocks(trajectory_lengths)
    elseif mode === :auto
        len = 0
        for c in axes(sums, 1)
            x = view(sums, c, :) ./ view(counts, c, :)
            l = optimal_block_length(x)
            l === nothing || (len = max(len, l))
        end
        if len == 0 || n ÷ len < 2
            @warn "automatic block-size selection is not possible ($n frames); using $n_blocks equal blocks" maxlog = 1
            return equal_blocks(n, n_blocks)
        end
        return blocks_of_length(n, len)
    else
        throw(ArgumentError("unknown blocking mode :$mode (use :equal, :trajectory or :auto)"))
    end
end
