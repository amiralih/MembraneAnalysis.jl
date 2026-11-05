# A small, dependency-free reader for the parts of a PDB file that the package needs.

"""
    PDBAtoms

Atom records read from a PDB file (structure-of-arrays layout). Atoms are numbered `1:N` in the
order in which they appear in the file; this order must match the order of atoms in the
trajectories.

Fields: `name`, `resname`, `chain`, `resnum` (residue number as written in the file),
`residue` (sequential residue counter starting at 1), `segname`, `x`, `y`, `z`, and `box`
(the `CRYST1` box lengths, or `nothing`).
"""
struct PDBAtoms
    name::Vector{String}
    resname::Vector{String}
    chain::Vector{String}
    resnum::Vector{Int}
    residue::Vector{Int}
    segname::Vector{String}
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    box::Union{Nothing,NTuple{3,Float64}}
end

Base.length(a::PDBAtoms) = length(a.name)

# columns are 1-based and inclusive, as in the PDB format description (ASCII lines only)
function _field(line::String, a::Int, b::Int)
    n = ncodeunits(line)
    n >= a || return ""
    return String(strip(SubString(line, a, min(b, n))))
end

"""
    read_pdb(path) -> PDBAtoms

Read `ATOM`/`HETATM` records (and the `CRYST1` box) of the first model in a PDB file.
Four-character residue names (columns 18–21), five-digit residue numbers (columns 23–27) and
segment names (columns 73–76), as written by CHARMM-GUI, are supported.
"""
function read_pdb(path::AbstractString)
    isfile(path) || throw(ArgumentError("PDB file not found: $path"))
    atoms = open(io -> _parse_pdb(io, String(path)), path)
    length(atoms) == 0 && throw(ArgumentError("no ATOM/HETATM records found in $path"))
    return atoms
end

function _parse_pdb(io::IO, path::String)
    name = String[]; resname = String[]; chain = String[]; segname = String[]
    resnum = Int[]; residue = Int[]
    x = Float64[]; y = Float64[]; z = Float64[]
    box = nothing
    previous_key = ("", "", "", "")
    residue_counter = 0
    line_number = 0
    for raw in eachline(io)
        line_number += 1
        if startswith(raw, "ENDMDL")
            break
        elseif startswith(raw, "CRYST1")
            line = String(raw)
            isascii(line) || continue
            lengths = (tryparse(Float64, _field(line, 7, 15)), tryparse(Float64, _field(line, 16, 24)),
                       tryparse(Float64, _field(line, 25, 33)))
            if all(!isnothing, lengths)
                box = (lengths[1], lengths[2], lengths[3])
            end
        elseif startswith(raw, "ATOM") || startswith(raw, "HETATM")
            line = String(raw)
            isascii(line) || throw(ArgumentError("$path, line $line_number: non-ASCII characters in an ATOM record"))
            ncodeunits(line) >= 54 ||
                throw(ArgumentError("$path, line $line_number: ATOM record is too short to contain coordinates"))
            xi = tryparse(Float64, _field(line, 31, 38))
            yi = tryparse(Float64, _field(line, 39, 46))
            zi = tryparse(Float64, _field(line, 47, 54))
            (xi === nothing || yi === nothing || zi === nothing) &&
                throw(ArgumentError("$path, line $line_number: could not read coordinates"))
            # column 27 is the insertion code in standard PDB files and the fifth digit of the
            # residue number in CHARMM-GUI files; both distinguish residues
            resnum_text = _field(line, 23, 27)
            key = (resnum_text, _field(line, 18, 21), _field(line, 22, 22), _field(line, 73, 76))
            if key != previous_key
                residue_counter += 1
                previous_key = key
            end
            push!(name, _field(line, 13, 16))
            push!(resname, key[2])
            push!(chain, key[3])
            push!(segname, key[4])
            push!(resnum, something(tryparse(Int, resnum_text), tryparse(Int, _field(line, 23, 26)), 0))
            push!(residue, residue_counter)
            push!(x, xi); push!(y, yi); push!(z, zi)
        end
    end
    return PDBAtoms(name, resname, chain, resnum, residue, segname, x, y, z, box)
end
