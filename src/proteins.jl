# Protein (peptide) definitions.

"""
    Protein(name; segments=String[], chains=String[], segment_prefix="PRO",
            backbone_atoms=["BB", "CA"], leaflet=:auto)

Description of a protein or peptide embedded in (or bound to) the membrane.

Protein atoms are selected from the PDB file by exact segment names (`segments`, PDB columns
73–76), by chain identifiers (`chains`), or — when neither is given — by a segment-name prefix
(`segment_prefix`, `"PRO"` by default, which matches CHARMM-GUI segment names `PROA`, `PROB`, …).
Several copies of the same protein (e.g. `PROA` and `PROB`) are treated together: per-residue
results are averaged over all copies.

Residues are identified by their residue number in the PDB file. Each residue is represented by
its backbone atom(s) (`backbone_atoms`; `BB` for Martini, `CA` for all-atom).

`leaflet` sets the sign convention used for curvature and heights of the protein:
`:upper` (+1), `:lower` (−1) or `:auto`, which uses the position of every backbone atom in the
PDB structure relative to the bilayer midplane of the structure.

# Example

```julia
peptide = Protein("PEP"; segments=["PROA"])
```
"""
struct Protein
    name::String
    segments::Vector{String}
    chains::Vector{String}
    segment_prefix::String
    backbone_atoms::Vector{String}
    leaflet::Symbol
end

function Protein(name::AbstractString; segments::AbstractVector=String[],
                 chains::AbstractVector=String[], segment_prefix::AbstractString="PRO",
                 backbone_atoms::AbstractVector=["BB", "CA"], leaflet::Symbol=:auto)
    leaflet in (:auto, :upper, :lower) ||
        throw(ArgumentError("Protein $name: leaflet must be :auto, :upper or :lower (got :$leaflet)"))
    isempty(backbone_atoms) && throw(ArgumentError("Protein $name: backbone_atoms is empty"))
    return Protein(String(name), String[String(s) for s in segments],
                   String[String(c) for c in chains], String(segment_prefix),
                   String[String(a) for a in backbone_atoms], leaflet)
end

function Base.show(io::IO, p::Protein)
    sel = !isempty(p.segments) ? "segments=$(p.segments)" :
          !isempty(p.chains) ? "chains=$(p.chains)" : "segment_prefix=\"$(p.segment_prefix)\""
    print(io, "Protein(\"", p.name, "\", ", sel, ", backbone_atoms=", p.backbone_atoms,
          ", leaflet=:", p.leaflet, ")")
end
