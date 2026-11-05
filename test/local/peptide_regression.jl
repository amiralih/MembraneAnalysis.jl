# Regression check of v0.4 against results of the experimental (MA_Hmin_pepH) code for the
# peptide-in-POPC sample (data_MA2026jun24). Not part of `Pkg.test()` because the data are not
# public. Run from the package directory with
#
#   julia --project=. test/local/peptide_regression.jl /path/to/data_MA2026jun24
#
# The data directory must contain step5_charmm2gmx.pdb, step7_1r.xtc and the reference outputs
# in t_1r/ (produced by as_newLs_r.jl, i.e. fluctuation_spectrum_hmin_allC with a 5 Å mask).
# Means must agree to within the Float32 precision of the old code; standard errors differ
# because v0.4 forms blocks of whole frames.

using MembraneAnalysis
using DelimitedFiles
using Test

data = length(ARGS) >= 1 ? ARGS[1] : get(ENV, "MA_PEPTIDE_DATA", "")
isdir(data) || error("usage: julia --project=. test/local/peptide_regression.jl /path/to/data_MA2026jun24")
reference = joinpath(data, "t_1r")
out = mktempdir()

params = joinpath(out, "parameters.toml")
write(params, """
    [[lipids]]
    name = "POPC"
    head_atom = "NC3"
    ref_atom = "C1A"
    tail_atom = "C4B"
    n_atoms = 12

    [[proteins]]
    name = "PEP"

    [surface]
    grid_spacing = 15.0
    mask_radius = 5.0

    [spectrum]
    q_min = 0.01
    q_max = 0.20
    n_bins = 19

    [[analysis]]
    type = "box"

    [[analysis]]
    type = "lipid_heights"

    [[analysis]]
    type = "protein_heights"
    residues = "1:8"

    [[analysis]]
    type = "lipid_curvature_spectrum"

    [[analysis]]
    type = "lipid_thickness_spectrum"

    [[analysis]]
    type = "protein_curvature_spectrum"
    residue = 4

    [[analysis]]
    type = "protein_curvature_spectrum"
    residue = 4
    direction = "x"

    [[analysis]]
    type = "protein_curvature_spectrum"
    residue = 4
    direction = "y"

    [[analysis]]
    type = "protein_thickness_spectrum"
    residue = 4

    [[analysis]]
    type = "protein_curvature_spectrum"
    residue = 5
    label = "cqs_fixTCB"

    [[analysis]]
    type = "protein_curvature_spectrum"
    residue = 5
    direction = "x"
    label = "cqs_X_fixTCB"

    [[analysis]]
    type = "protein_curvature_spectrum"
    residue = 5
    direction = "y"
    label = "cqs_Y_fixTCB"
    """)

run_analysis(params, joinpath(data, "step5_charmm2gmx.pdb"), [joinpath(data, "step7_1r.xtc")];
             output_dir=out, nprocs=1)

function compare(name; rtol=1e-3)
    old = readdlm(joinpath(reference, name))
    new = readdlm(joinpath(out, name); comments=true)
    a = Float64.(new[:, 2])
    b = Float64.(old[:, 2])
    ok = isfinite.(b)
    scale = maximum(abs, b[ok])
    dev = maximum(abs.(a[ok] .- b[ok]))
    println(rpad(name, 24), "  max |Δmean| = ", round(dev; sigdigits=3), "  (max |mean| = ",
            round(scale; sigdigits=3), ")")
    # Bins whose modes move in and out of the bin as the box fluctuates were NaN in the old code
    # (any frame without modes made the whole bin NaN); v0.4 averages them over the frames in
    # which they contain modes. Only bins with a value in the old output are compared.
    extra = findall(isfinite.(a) .& .!ok)
    isempty(extra) || println("    bins populated in only some frames (NaN in the old code): q = ",
                              join(new[extra, 1], ", "))
    @test all(isfinite, a[ok])
    @test all(isapprox.(a[ok], b[ok]; rtol=rtol, atol=rtol * scale))
end

@testset "peptide sample (data_MA2026jun24)" begin
    for f in ["POPC_zs.dat", "PEP_zs.dat", "POPC_cqs.dat", "POPC_tqs.dat", "PEP_cqs.dat",
              "PEP_cqs_X.dat", "PEP_cqs_Y.dat", "PEP_tqs.dat", "PEP_cqs_fixTCB.dat",
              "PEP_cqs_X_fixTCB.dat", "PEP_cqs_Y_fixTCB.dat"]
        compare(f)
    end
    @test readdlm(joinpath(out, "KA.dat"); comments=true)[1] ≈ readdlm(joinpath(reference, "KA.dat"))[1] rtol = 1e-8
    @test readdlm(joinpath(out, "A.dat"); comments=true)[1] ≈ readdlm(joinpath(reference, "A.dat"))[1] rtol = 1e-4
end
println("outputs of this check: ", out)
