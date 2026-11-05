@testset "parameters file" begin
    dir = mktempdir()
    function cfg_from(text)
        path = joinpath(dir, "params_$(hash(text)).toml")
        write(path, text)
        return load_parameters(path)
    end

    base = """
    title = "test"

    [[lipids]]
    preset = "POPC_m2"
    ref_atom = "C1A"

    [[lipids]]
    name = "CHOL"
    head_atom = "ROH"
    ref_atom = "R1"
    tail_atom = "C2"
    n_atoms = 8

    [[proteins]]
    name = "PEP"
    segments = ["PROA"]

    [surface]
    grid_spacing = 15
    mask_radius = 5.0

    [statistics]
    blocking = "auto"

    [[analysis]]
    type = "box"

    [[analysis]]
    type = "protein_curvature_spectrum"
    residue = 4
    direction = "x"
    q_max = 0.2
    n_bins = 19
    """
    cfg = cfg_from(base)
    @test cfg.title == "test"
    @test length(cfg.lipids) == 2
    @test cfg.lipids[1].name == "POPC" && cfg.lipids[1].ref_atom == "C1A" && cfg.lipids[1].n_atoms == 12
    @test cfg.lipids[2].name == "CHOL" && cfg.lipids[2].n_atoms == 8
    @test cfg.proteins[1].segments == ["PROA"]
    @test cfg.grid_spacing == 15.0
    @test cfg.fill === :min_curvature          # default when masking
    @test cfg.mask_proteins == ["PEP"]         # all proteins by default
    @test cfg.blocking === :auto
    @test cfg.n_blocks == 10
    @test cfg.save_spectra
    @test length(cfg.analyses) == 2
    @test cfg.spectrum_defaults["q_min"] == 0.01

    # every TOML value that is written back must be valid TOML
    io = IOBuffer()
    MembraneAnalysis.TOML.print(io, MembraneAnalysis.config_dict(cfg); sorted=true)
    @test MembraneAnalysis.TOML.parse(String(take!(io)))["lipids"][1]["name"] == "POPC"

    # residue lists
    @test MembraneAnalysis._parse_int_list("1:4") == [1, 2, 3, 4]
    @test MembraneAnalysis._parse_int_list("1-3, 7") == [1, 2, 3, 7]
    @test MembraneAnalysis._parse_int_list([2, 5]) == [2, 5]

    # errors
    @test_throws ParameterError cfg_from(replace(base, "title = \"test\"" => "titel = \"test\""))
    @test_throws ParameterError cfg_from(replace(base, "direction = \"x\"" => "direction = \"z\""))
    @test_throws ParameterError cfg_from(replace(base, "type = \"box\"" => "type = \"boxes\""))
    @test_throws ParameterError cfg_from(replace(base, "preset = \"POPC_m2\"" => "preset = \"POPC_m9\""))
    @test_throws ParameterError cfg_from(replace(base, "blocking = \"auto\"" => "blocking = \"weekly\""))
    @test_throws ParameterError cfg_from(replace(base, "mask_radius = 5.0" => "mask_radius = 5.0\nfill = \"neighbors\""))
    @test_throws ParameterError cfg_from(replace(base, "residue = 4" => "residue = 4\nresidues = [1, 2]"))
    @test_throws ParameterError cfg_from(replace(base, "n_atoms = 8" => "n_atoms = \"eight\""))
    @test_throws ParameterError cfg_from(replace(base, "residue = 4\n" => ""))
    @test_throws ParameterError cfg_from(replace(base, "type = \"box\"" => "type = \"undulation\""))
    @test_throws ParameterError cfg_from(replace(base, "type = \"box\"" => "type = \"undulation\"\nq_max = 0.1\nmodel = \"XY\""))
    @test cfg_from(replace(base, "type = \"box\"" => "type = \"undulation\"\nq_max = 0.1")) isa AnalysisConfig
    @test_throws ParameterError cfg_from("""
        [[lipids]]
        preset = "DOPC_m2"
        """)
    @test_throws ParameterError load_parameters(joinpath(dir, "does_not_exist.toml"))
    bad = joinpath(dir, "bad.toml")
    write(bad, "this is = = not toml")
    @test_throws ParameterError load_parameters(bad)

    @test "tcb" in MembraneAnalysis.analysis_types()

    # command-line parsing
    o = MembraneAnalysis.parse_cli(["p.toml", "s.pdb", "t1.xtc", "t2.xtc", "-o", "out", "-n", "4", "--resume"])
    @test o.parameters == "p.toml" && o.structure == "s.pdb" && o.trajectories == ["t1.xtc", "t2.xtc"]
    @test o.output == "out" && o.nprocs == 4 && o.resume && !o.check
    o = MembraneAnalysis.parse_cli(["--nprocs=2", "--check", "p.toml", "s.pdb", "t.xtc"])
    @test o.nprocs == 2 && o.check
    @test MembraneAnalysis.parse_cli(["--help"]).action === :help
    @test_throws ParameterError MembraneAnalysis.parse_cli(["p.toml", "s.pdb"])
    @test_throws ParameterError MembraneAnalysis.parse_cli(["p.toml", "s.pdb", "t.xtc", "--bogus"])
    @test_throws ParameterError MembraneAnalysis.parse_cli(["p.toml", "s.pdb", "t.xtc", "-n", "zero"])
    @test MembraneAnalysis.cli_main(["--version"]) == 0
    @test MembraneAnalysis.cli_main(["p.toml"]) == 1
end
