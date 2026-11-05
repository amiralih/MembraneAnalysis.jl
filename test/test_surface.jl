@testset "surface reconstruction" begin
    @test MembraneAnalysis.fft_index.(1:6, 6) == [0, 1, 2, -3, -2, -1]
    @test MembraneAnalysis.fft_index.(1:5, 5) == [0, 1, 2, -2, -1]

    w = MembraneAnalysis.bending_weights(8, 6)
    @test w[2, 1] == 1 && w[8, 1] == 1 && w[1, 2] == 1 && w[2, 2] == 4
    @test w[5, 1] == 0          # Nyquist mode excluded
    @test all(w[i, j] == w[mod(-(i - 1), 8) + 1, mod(-(j - 1), 6) + 1] for i in 1:8, j in 1:6)

    # neighbour filling: cells are filled in order, using already filled neighbours
    zsum = Float64[1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16]
    ns = ones(Int, 4, 4)
    for c in [(1, 1), (2, 1), (3, 3)]
        zsum[c...] = 0
        ns[c...] = 0
    end
    h = zeros(4, 4)
    MembraneAnalysis.fill_neighbors!(h, zsum, ns)
    @test h ≈ [19/3 2 3 4; 7 6 7 8; 9 10 11 12; 13 14 15 16]

    # curvature-minimising filling equals the mode-by-mode normal equations of the
    # experimental fluctuation_spectrum_hmin functions
    function sum_q4s(H)
        s = 0.0im
        nx, ny = size(H)
        for i in -(nx ÷ 2 - 1):(nx ÷ 2 - 1), j in -(ny ÷ 2 - 1):(ny ÷ 2 - 1)
            s += H[mod1(i + 1, nx), mod1(j + 1, ny)] * (i^2 + j^2)^2
        end
        return s
    end
    function bruteforce_fill(h, cells)
        h0 = copy(h)
        for c in cells
            h0[c] = 0.0
        end
        hf0 = MembraneAnalysis.fft(h0)
        E = map(cells) do c
            e = zeros(size(h))
            e[c] = 1.0
            MembraneAnalysis.fft(e)
        end
        M = length(cells)
        X = [real(sum_q4s(E[i] .* conj.(E[j]))) for i in 1:M, j in 1:M]
        Y = [-real(sum_q4s(hf0 .* conj.(E[i]))) for i in 1:M]
        out = copy(h)
        out[cells] .= X \ Y
        return out
    end
    for (nx, ny, nmask) in [(14, 14, 20), (13, 15, 30), (10, 12, 8)]
        h = reshape(lcg_uniform(nx * ny; seed=nx * ny) .* 10 .- 5, nx, ny)
        mask = falses(nx, ny)
        mask[sort!(unique(round.(Int, lcg_uniform(nmask; seed=7) .* (nx * ny - 1))) .+ 1)] .= true
        ws = SurfaceWorkspace(SurfaceSettings(grid_size=(nx, ny)), 0)
        a = MembraneAnalysis.fill_min_curvature!(copy(h), mask, ws)
        b = bruteforce_fill(h, findall(mask))
        @test maximum(abs, a .- b) < 1e-8
        @test a[.!mask] == h[.!mask]
        # the result is a minimum: perturbing a masked cell increases the bending energy
        energy(x) = sum(ws.q4 .* abs2.(MembraneAnalysis.fft(x)))
        c = findfirst(mask)
        for δ in (-1e-3, 1e-3)
            p = copy(a)
            p[c] += δ
            @test energy(p) > energy(a)
        end
    end
    ws = SurfaceWorkspace(SurfaceSettings(grid_size=(6, 6)), 0)
    h = reshape(collect(1.0:36.0), 6, 6)
    @test MembraneAnalysis.fill_min_curvature!(copy(h), falses(6, 6), ws) == h

    @test_throws ArgumentError SurfaceSettings(grid_size=(10, 10), fill=:neighbors, mask_radius=5.0)
    @test_throws ArgumentError SurfaceSettings(grid_size=(2, 10))
    @test MembraneAnalysis.grid_size_for(198.66, 198.66, 15.0) == (14, 14)
end
