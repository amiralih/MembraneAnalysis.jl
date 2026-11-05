using MembraneAnalysis
using Test
using DelimitedFiles

const SAMPLE_DIR = joinpath(@__DIR__, "sample_files")

# deterministic pseudo-random numbers (independent of the Julia version's RNG)
function lcg_uniform(n::Integer; seed::Integer=12345)
    x = UInt64(seed)
    out = zeros(n)
    for i in 1:n
        x = x * 0x5851f42d4c957f2d + 0x14057b7ef767814f
        out[i] = (x >> 11) / 2.0^53
    end
    return out
end

@testset "MembraneAnalysis.jl" begin
    include("test_statistics.jl")
    include("test_fourier.jl")
    include("test_surface.jl")
    include("test_system.jl")
    include("test_parameters.jl")
    include("test_pipeline.jl")
end
