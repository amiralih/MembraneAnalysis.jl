# Lipid species definitions and presets.

"""
    Lipid(name, head_atom, ref_atom, tail_atom, n_atoms; height_atoms=String[])

Description of a lipid species present in the membrane.

# Fields

- `name`: residue name of the lipid in the PDB file (e.g. `"POPC"`).
- `head_atom`: name of the first atom of every lipid molecule in the PDB file. The leaflet of a
  lipid is decided from the height of this atom.
- `ref_atom`: default reference atom. Its position defines the membrane surface (fluctuation
  spectrum) and the position of the lipid in the lipid spectra analyses.
- `tail_atom`: an atom at the end of a lipid tail; the mean height of all tail atoms is the
  bilayer midplane.
- `n_atoms`: number of atoms in one lipid molecule. Atoms of one molecule must be contiguous in
  the PDB file, starting with `head_atom`.
- `height_atoms`: atoms used by the atom-height, sampled-curvature and TCB analyses. When empty
  (default) all atoms whose name starts with `C`, `D`, `G` or `X` are used (carbon atoms in
  all-atom force fields, glycerol/tail beads in Martini).

# Example

```julia
POPC = Lipid("POPC", "NC3", "C1A", "C4B", 12)
```
"""
struct Lipid
    name::String
    head_atom::String
    ref_atom::String
    tail_atom::String
    n_atoms::Int
    height_atoms::Vector{String}

    function Lipid(name::AbstractString, head_atom::AbstractString, ref_atom::AbstractString,
                   tail_atom::AbstractString, n_atoms::Integer, height_atoms::AbstractVector)
        n_atoms > 0 || throw(ArgumentError("Lipid $name: n_atoms must be positive (got $n_atoms)"))
        return new(String(name), String(head_atom), String(ref_atom), String(tail_atom),
                   Int(n_atoms), String[String(a) for a in height_atoms])
    end
end

function Lipid(name::AbstractString, head_atom::AbstractString, ref_atom::AbstractString,
               tail_atom::AbstractString, n_atoms::Integer; height_atoms::AbstractVector=String[])
    return Lipid(name, head_atom, ref_atom, tail_atom, n_atoms, height_atoms)
end

"""
    with_ref_atom(lipid, atom)

Return a copy of `lipid` that uses `atom` as its reference atom.
"""
with_ref_atom(l::Lipid, atom::AbstractString) =
    Lipid(l.name, l.head_atom, atom, l.tail_atom, l.n_atoms, l.height_atoms)

function Base.show(io::IO, l::Lipid)
    print(io, "Lipid(\"", l.name, "\", head=\"", l.head_atom, "\", ref=\"", l.ref_atom,
          "\", tail=\"", l.tail_atom, "\", n_atoms=", l.n_atoms, ")")
end

# --------------------------------------------------------------------------- presets
# All-atom (CHARMM36)

"All-atom (CHARMM36) POPC."
const POPC_aa = Lipid("POPC", "N", "C37", "C316", 134)
"All-atom (CHARMM36) cholesterol."
const Chol_aa = Lipid("CHL1", "C3", "C11", "C27", 74)
"All-atom (CHARMM36) DMPC."
const DMPC_aa = Lipid("DMPC", "N", "C3", "C314", 118)
"All-atom (CHARMM36) DOPC."
const DOPC_aa = Lipid("DOPC", "N", "C24", "C318", 138)
"All-atom (CHARMM36) SOPC."
const SOPC_aa = Lipid("SOPC", "N", "C22", "C318", 140)

# Martini 2

"Martini 2 POPC."
const POPC_m2 = Lipid("POPC", "NC3", "C1B", "C4B", 12)
"Martini 2 DOPE."
const DOPE_m2 = Lipid("DOPE", "NH3", "D2A", "C4B", 12)
"Martini 2 DOPC."
const DOPC_m2 = Lipid("DOPC", "NC3", "C1B", "C4B", 12)

"""
Lipid presets that can be referred to by name in a parameters file (`preset = "POPC_m2"`).
"""
const LIPID_PRESETS = Dict{String,Lipid}(
    "POPC_aa" => POPC_aa,
    "Chol_aa" => Chol_aa,
    "DMPC_aa" => DMPC_aa,
    "DOPC_aa" => DOPC_aa,
    "SOPC_aa" => SOPC_aa,
    "POPC_m2" => POPC_m2,
    "DOPE_m2" => DOPE_m2,
    "DOPC_m2" => DOPC_m2,
)
