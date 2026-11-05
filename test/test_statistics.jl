@testset "statistics" begin
    x = sin.(1:1000)

    # value from the test suite of earlier versions
    @test blocking_error(x) ≈ 0.0012181125344636647

    @test equal_blocks(25, 10) == [(2b - 1):(2b) for b in 1:10]
    @test equal_blocks(100, 4) == [1:25, 26:50, 51:75, 76:100]
    @test trajectory_blocks([3, 2, 4]) == [1:3, 4:5, 6:9]
    @test MembraneAnalysis.blocks_of_length(10, 4) == [1:4, 5:8]

    # one component with unit counts: block statistics reduce to blocking_error
    bs = BlockStats(["x"], reshape(x, 1, :), ones(1, 1000), equal_blocks(1000, 10))
    @test block_means(bs)[1] ≈ sum(x) / 1000
    @test block_sems(bs)[1] ≈ blocking_error(x)

    # for equal blocks covering all frames, the jackknife error of the mean is the block SEM
    est, err = jackknife(v -> v[1], bs)
    @test est ≈ sum(x) / 1000
    @test err ≈ blocking_error(x)

    # vector-valued estimators
    est2, err2 = jackknife(v -> [v[1], 2v[1]], bs)
    @test est2 ≈ [1, 2] .* (sum(x) / 1000)
    @test err2 ≈ [1, 2] .* blocking_error(x)

    # ratio estimator: components with varying counts, never-sampled component gives NaN
    sums = [1.0 2.0 3.0 4.0; 0.0 0.0 0.0 0.0]
    counts = [1.0 1.0 2.0 2.0; 0.0 0.0 0.0 0.0]
    bs2 = BlockStats(["a", "b"], sums, counts, [1:2, 3:4])
    m = block_means(bs2)
    @test m[1] ≈ 10 / 6
    @test isnan(m[2])
    @test block_sems(bs2)[1] ≈ 0.125     # block means 1.5 and 1.75
    @test isnan(block_sems(bs2)[2])

    # automatic block size: white noise needs no blocking, an AR(1) process does
    u = lcg_uniform(4096) .- 0.5
    @test optimal_block_length(u) == 1
    ar = zeros(4096)
    for t in 2:4096
        ar[t] = 0.9 * ar[t - 1] + u[t]
    end
    @test optimal_block_length(ar) == 32
    @test optimal_block_length(u[1:10]) === nothing
    @test optimal_block_length(ones(100)) === nothing

    @test MembraneAnalysis.choose_blocks(:equal, zeros(2, 50), ones(2, 50), [50]; n_blocks=5) ==
          equal_blocks(50, 5)
    @test MembraneAnalysis.choose_blocks(:trajectory, zeros(2, 50), ones(2, 50), [20, 30]) ==
          [1:20, 21:50]
    blocks = MembraneAnalysis.choose_blocks(:auto, reshape(ar, 1, :), ones(1, 4096), [4096])
    @test length(first(blocks)) == 32
end
