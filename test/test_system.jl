@testset "structure and system" begin
    # a tiny synthetic system: a 3-residue peptide and 4 three-atom lipids
    function pdb_line(serial, name, resname, resnum, x, y, z, seg)
        f(v) = lpad(string(round(v; digits=3)), 8)
        return "ATOM  " * lpad(serial, 5) * " " * rpad(name, 4) * " " * rpad(resname, 4) * " " *
               lpad(resnum, 4) * "    " * f(x) * f(y) * f(z) * "  1.00  0.00      " * rpad(seg, 4)
    end
    lines = ["CRYST1   50.000   40.000   80.000  90.00  90.00  90.00 P 1           1"]
    serial = 0
    for r in 1:3, name in ("BB", "SC1")
        serial += 1
        push!(lines, pdb_line(serial, name, "LEU", r, 10.0 + 3r, 20.0, 55.0, "PROA"))
    end
    for (m, (zhead, ztail)) in enumerate([(60.0, 45.0), (61.0, 44.0), (20.0, 35.0), (21.0, 36.0)])
        for (name, z) in (("H1", zhead), ("C1", (zhead + ztail) / 2), ("C2", ztail))
            serial += 1
            push!(lines, pdb_line(serial, name, "LIPA", m, 5.0 * m, 7.0, z, "MEMB"))
        end
    end
    push!(lines, "END")
    path = joinpath(mktempdir(), "mini.pdb")
    write(path, join(lines, "\n") * "\n")

    atoms = read_pdb(path)
    @test length(atoms) == 18
    @test atoms.box == (50.0, 40.0, 80.0)
    @test atoms.name[1:2] == ["BB", "SC1"]
    @test atoms.segname[1] == "PROA" && atoms.segname[end] == "MEMB"
    @test atoms.resname[7] == "LIPA" && atoms.resnum[7] == 1
    @test atoms.residue[end] == 7
    @test atoms.z[7] ≈ 60.0

    lipa = Lipid("LIPA", "H1", "C1", "C2", 3)
    pep = Protein("PEP"; segments=["PROA"])
    s = MembraneSystem(path, [lipa], [pep])
    @test s.n_atoms == 18
    @test s.mol_head == [7, 10, 13, 16]
    @test s.static_leaflet == Int8[1, 1, -1, -1]
    @test s.tail_atoms == [9, 12, 15, 18]
    @test s.pdb_midplane ≈ 40.0
    @test s.atom_mol[7:9] == [1, 1, 1] && s.atom_mol[1] == 0
    @test s.surface_atoms == [8, 11, 14, 17]
    @test MembraneAnalysis.protein_residues(s, 1) == [1, 2, 3]
    g = MembraneAnalysis.protein_residue_group(s, 1, 2)
    @test g.atoms == [3] && g.signs == Int8[1]
    @test [x.label for x in MembraneAnalysis.height_atom_groups(s, 1)] == ["C1", "C2"]
    s_lower = MembraneSystem(path, [lipa], [Protein("PEP"; segments=["PROA"], leaflet=:lower)])
    @test MembraneAnalysis.protein_residue_group(s_lower, 1, 2).signs == Int8[-1]
    @test sprint(show, s) isa String
    @test sprint(show, MIME("text/plain"), s) isa String

    # default protein selection by segment-name prefix "PRO"
    @test MembraneSystem(path, [lipa], [Protein("PEP")]).protein_backbone[1] == [1, 3, 5]

    @test_throws ArgumentError MembraneSystem(path, [Lipid("LIPA", "H1", "C1", "C2", 4)])
    @test_throws ArgumentError MembraneSystem(path, [Lipid("LIPA", "C1", "C1", "C2", 2)])
    @test_throws ArgumentError MembraneSystem(path, [Lipid("POPC", "NC3", "C1A", "C4B", 12)])
    @test_throws ArgumentError MembraneSystem(path, [Lipid("LIPA", "H1", "XX", "C2", 3)])
    @test_throws ArgumentError MembraneSystem(path, [lipa], [Protein("P2"; segments=["PROZ"])])
    @test_throws ArgumentError MembraneAnalysis.protein_residue_group(s, 1, 9)

    # five-digit residue numbers (columns 23-27) as written by CHARMM-GUI
    waters = joinpath(mktempdir(), "waters.pdb")
    write(waters, """
        ATOM  33999  W   W    9999     216.930 219.569  78.777  1.00  0.00      WAT
        ATOM  34000  W   W    10000    186.290 248.179  73.187  1.00  0.00      WAT
        ATOM  34001  W   W    10001    197.600 243.009  75.667  1.00  0.00      WAT
        END
        """)
    w = read_pdb(waters)
    @test w.resnum == [9999, 10000, 10001]
    @test w.residue == [1, 2, 3]
    @test w.z ≈ [78.777, 73.187, 75.667]

    # the public sample system (Martini DOPE/DOPC)
    sample = MembraneSystem(joinpath(SAMPLE_DIR, "s_1.pdb"), [DOPE_m2, DOPC_m2])
    @test sample.n_atoms == 45509
    @test length(sample.mol_head) == 2000
    @test length(sample.tail_atoms) == 2000
    @test length(sample.ref_atoms[1]) == 200 && length(sample.ref_atoms[2]) == 1800
end
