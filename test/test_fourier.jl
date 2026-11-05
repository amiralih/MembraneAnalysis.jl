@testset "Fourier evaluation" begin
    # index pairs (test of earlier versions)
    values, pairs = MembraneAnalysis.get_index_pairs(2)
    @test values ≈ [0.0, 1.0, 1.4142135623730951, 2.0]
    @test pairs == [[(0, 0)], [(0, 1), (1, 0)], [(1, 1)], [(0, 2), (2, 0)]]

    b = spectrum_bins(0.01, 0.21, 10)
    @test length(b) == 10
    @test b.lo[1] == 0.01
    @test b.hi[end] ≈ 0.21
    @test b.centers ≈ 0.02:0.02:0.20
    @test_throws ArgumentError spectrum_bins(0.2, 0.1, 10)

    # a single cosine mode h(x, y) = a cos(qx x): with h_q normalised so that
    # h(r) = (1 / (Lx Ly)) ∑ h_q e^{iq·r}, the curvature per mode is qx² a cos(qx x)
    nx, ny = 12, 10
    Lx, Ly = 180.0, 150.0
    a = 2.5
    qx = 2π / Lx
    hq = zeros(ComplexF64, nx, ny)
    hq[2, 1] = a * Lx * Ly / 2       # mode (1, 0)
    hq[nx, 1] = a * Lx * Ly / 2      # mode (-1, 0)
    # q_max = 1.1 qx keeps only the modes (±1, 0) (2π/Ly > 1.1 qx); for direction = :x the
    # modes (0, ±1) are excluded anyway
    for (x, y) in [(0.0, 0.0), (37.0, 12.0), (95.5, 140.0)]
        @test curvature_at(hq, Lx, Ly, x, y; q_max=1.1qx) ≈ qx^2 * a * cos(qx * x)
        @test curvature_at(hq, Lx, Ly, x, y; q_max=1.5qx, direction=:x) ≈ qx^2 * a * cos(qx * x)
        @test thickness_at(hq, Lx, Ly, x, y; q_max=1.1qx) ≈ a * cos(qx * x)
        # with the (empty) modes (0, ±1) included the value per mode is halved
        @test curvature_at(hq, Lx, Ly, x, y; q_max=1.5qx) ≈ qx^2 * a * cos(qx * x) / 2
    end
    # modes along y exist in the shell but carry no amplitude
    @test curvature_at(hq, Lx, Ly, 10.0, 10.0; q_max=1.5qx, direction=:y) ≈ 0.0 atol = 1e-12
    # a shell without modes gives NaN
    @test isnan(curvature_at(hq, Lx, Ly, 10.0, 10.0; q_min=0.0, q_max=0.5qx))

    # bins: modes are assigned to at most one bin and counted
    f = FourierField(hq, Lx, Ly, spectrum_bins(0.0, 0.2, 8), :all, :curvature)
    @test length(f.coef) == sum(f.n_modes)
    @test all(1 .<= f.bin .<= 8)
    fx = FourierField(hq, Lx, Ly, spectrum_bins(0.0, 0.2, 8), :x, :curvature)
    @test all(==(0), fx.j)
    fy = FourierField(hq, Lx, Ly, spectrum_bins(0.0, 0.2, 8), :y, :curvature)
    @test all(==(0), fy.i)
end
