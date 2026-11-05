@testset "pipeline (regression against v0.3.2 results)" begin
    pdb = joinpath(SAMPLE_DIR, "s_1.pdb")
    traj = joinpath(SAMPLE_DIR, "t_1.xtc")
    ref(name) = readdlm(joinpath(SAMPLE_DIR, name))
    new(dir, name) = readdlm(joinpath(dir, name); comments=true)

    dir = mktempdir()
    params = joinpath(dir, "parameters.toml")
    write(params, """
        [[lipids]]
        preset = "DOPE_m2"

        [[lipids]]
        preset = "DOPC_m2"

        [surface]
        grid_spacing = 15.0

        [[analysis]]
        type = "box"

        [[analysis]]
        type = "undulation"
        q_max = 0.08

        [[analysis]]
        type = "lipid_heights"

        [[analysis]]
        type = "lipid_sampled_curvature"
        q_max = 0.08

        [[analysis]]
        type = "lipid_curvature_spectrum"
        q_min = 0.01
        q_max = 0.21
        n_bins = 10

        [[analysis]]
        type = "lipid_thickness_spectrum"
        q_min = 0.01
        q_max = 0.21
        n_bins = 10

        [[analysis]]
        type = "lipid_curvature_spectrum"
        direction = "x"
        q_min = 0.01
        q_max = 0.21
        n_bins = 10

        [[analysis]]
        type = "tcb"
        q_max = 0.08
        z_cutoff = 20.0
        """)

    out1 = joinpath(dir, "serial")
    summary = run_analysis(params, pdb, [traj]; output_dir=out1, nprocs=1)

    # box, area and area expansion modulus
    @test Float64.(new(out1, "box_dims.dat")[:, 2]) ≈ vec(ref("box_dims.dat")) rtol = 1e-10
    @test new(out1, "A.dat")[1] ≈ ref("A.dat")[1] rtol = 1e-4      # now <Lx Ly> instead of <Lx><Ly>
    @test new(out1, "KA.dat")[1] ≈ ref("KA.dat")[1] rtol = 1e-6
    @test summary["box"]["KA"] ≈ ref("KA.dat")[1] rtol = 1e-6
    @test summary["box"]["KA_error"] > 0

    # undulation spectrum and bending modulus
    hq2 = new(out1, "hq2.dat")
    @test hq2[:, 1] ≈ ref("hq2.dat")[:, 1] rtol = 1e-10
    @test hq2[:, 2] ≈ ref("hq2.dat")[:, 2] rtol = 1e-4
    @test new(out1, "kc_hq2.dat")[1] ≈ ref("kc_hq2.dat")[1] rtol = 1e-3
    @test isfile(joinpath(out1, "tq2.dat"))

    for lipid in ("DOPE", "DOPC")
        zs = new(out1, "$(lipid)_zs.dat")
        @test zs[:, 1] == ref("$(lipid)_zs.dat")[:, 1]
        @test Float64.(zs[:, 2]) ≈ Float64.(ref("$(lipid)_zs.dat")[:, 2]) rtol = 1e-8
        cs = new(out1, "$(lipid)_cs.dat")
        @test Float64.(cs[:, 2]) ≈ Float64.(ref("$(lipid)_cs.dat")[:, 2]) rtol = 1e-4
        cqs = new(out1, "$(lipid)_cqs.dat")
        r = ref("$(lipid)_cqs.dat")
        @test cqs[:, 1] ≈ r[:, 1]
        ok = isfinite.(r[:, 2])
        @test isnan.(cqs[:, 2]) == .!ok
        @test cqs[ok, 2] ≈ r[ok, 2] rtol = 1e-4
        @test all(isfinite, new(out1, "$(lipid)_tqs.dat")[ok, 2])
        @test isfile(joinpath(out1, "$(lipid)_cqs_X.dat"))
        @test isfile(joinpath(out1, "$(lipid)_TCB_profile.dat"))
    end
    tcb = new(out1, "TCB.dat")
    @test tcb[1, 1] == "kc"
    @test Float64(tcb[1, 2]) ≈ ref("kc.dat")[1] rtol = 1e-4
    @test summary["tcb"]["kc"] ≈ ref("kc.dat")[1] rtol = 1e-4

    # per-trajectory time series
    ts = joinpath(out1, "timeseries", "001_t_1.h5")
    @test isfile(ts)
    @test isfile(joinpath(out1, "summary.toml"))
    @test isfile(joinpath(out1, "parameters_used.toml"))

    # resuming reuses the time series and gives identical results
    before = read(joinpath(out1, "DOPC_cqs.dat"), String)
    run_analysis(params, pdb, [traj]; output_dir=out1, nprocs=1, resume=true)
    @test read(joinpath(out1, "DOPC_cqs.dat"), String) == before

    # parallel processing with the trajectory split into chunks gives identical results
    out2 = joinpath(dir, "parallel")
    run_analysis(params, pdb, [traj]; output_dir=out2, nprocs=2, chunk_frames=4)
    for f in ("A.dat", "KA.dat", "hq2.dat", "DOPE_zs.dat", "DOPC_cs.dat", "DOPC_cqs.dat",
              "DOPE_tqs.dat", "TCB.dat")
        a = new(out1, f)
        b = new(out2, f)
        @test size(a) == size(b)
        @test all(x -> x[1] == x[2] || (x[1] isa Number && x[2] isa Number &&
                                        (isnan(x[1]) && isnan(x[2]) || isapprox(x[1], x[2]; rtol=1e-12))),
                  zip(a, b))
    end

    # checking the inputs only
    @test MembraneAnalysis.cli_main(["--check", "-n", "1", params, pdb, traj]) == 0

    # analyses that would overwrite each other's files are rejected
    clash = joinpath(dir, "clash.toml")
    write(clash, """
        [[lipids]]
        preset = "DOPC_m2"
        [[lipids]]
        preset = "DOPE_m2"
        [[analysis]]
        type = "lipid_heights"
        label = "cs"
        [[analysis]]
        type = "lipid_sampled_curvature"
        q_max = 0.08
        """)
    @test_throws ParameterError run_analysis(clash, pdb, traj; nprocs=1, check=true)
    # the same analysis type for different lipids needs no label
    subsets = joinpath(dir, "subsets.toml")
    write(subsets, """
        [[lipids]]
        preset = "DOPC_m2"
        [[lipids]]
        preset = "DOPE_m2"
        [[analysis]]
        type = "lipid_curvature_spectrum"
        lipids = ["DOPC"]
        [[analysis]]
        type = "lipid_curvature_spectrum"
        lipids = ["DOPE"]
        n_bins = 5
        """)
    @test run_analysis(subsets, pdb, traj; nprocs=1, check=true) isa Dict
end
