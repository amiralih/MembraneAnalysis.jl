# Writing result files.

"""
    OutputContext

Where and how results are written: the output `directory` and header lines (`header`)
describing the run, which are written at the top of every result file.
"""
struct OutputContext
    directory::String
    header::Vector{String}
end

_fmt(x::AbstractFloat) = isnan(x) ? "NaN" : string(Float64(x))
_fmt(x::Integer) = string(x)
_fmt(x) = string(x)

"""
    write_table(out, filename, description, column_names, columns...)

Write a tab-separated text file with `#` header lines (run information, `description` and the
column names) followed by one row per entry of the `columns`. Read it back with
`readdlm(file; comments=true)`.
"""
function write_table(out::OutputContext, filename::AbstractString, description::AbstractString,
                     column_names::Vector{String}, columns::AbstractVector...)
    n = length(first(columns))
    all(c -> length(c) == n, columns) || throw(DimensionMismatch("columns differ in length"))
    path = joinpath(out.directory, filename)
    open(path, "w") do io
        for line in out.header
            println(io, "# ", line)
        end
        println(io, "# ", description)
        println(io, "# ", join(column_names, '\t'))
        for r in 1:n
            println(io, join((_fmt(c[r]) for c in columns), '\t'))
        end
    end
    return path
end

"""
    stats_header(bs) -> String

One-line description of the blocks used for the error estimates.
"""
function stats_header(bs::BlockStats)
    lens = length.(bs.blocks)
    size_text = isempty(lens) ? "-" : minimum(lens) == maximum(lens) ? string(first(lens)) :
                "$(minimum(lens))–$(maximum(lens))"
    return "errors: standard error of the mean from $(n_blocks(bs)) blocks of $(size_text) " *
           "frames (block jackknife for derived quantities)"
end
