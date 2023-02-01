"""
    Lipid(
        name,
        head_atom,
        ref_atom,
        tail_atom,
        n_atoms)

A struct describing lipid species present in the system.

### Fields

* `name`: lipid name;
* `head_atom`: the head atom of the lipid (has to be the first atom in the PDB file);
* `ref_atom`: default reference atom for determining the position of the lipid;
* `tail_atom`: a tail atom of the lipid (does not have to be the last atom in the PDB file);
* `n_atoms`: number of atoms in the lipid.

"""
struct Lipid
    name
    head_atom
    ref_atom
    tail_atom
    n_atoms
end

# all-atom

POPC_aa = Lipid(
    "POPC",
    "N",
    "C37",
    "C316",
    134
)

Chol_aa = Lipid(
    "CHL1",
    "C3",
    "C11",
    "C27",
    74
)

DMPC_aa = Lipid(
    "DMPC",
    "N",
    "C3",
    "C314",
    118
)

DOPC_aa = Lipid(
    "DOPC",
    "N",
    "C24",
    "C318",
    138
)

SOPC_aa = Lipid(
    "SOPC",
    "N",
    "C22",
    "C318",
    140
)

# Martini 2

POPC_m2 = Lipid(
    "POPC",
    "NC3",
    "C1B",
    "C4B",
    12
)

DOPE_m2 = Lipid(
    "DOPE",
    "NH3",
    "D2A",
    "C4B",
    12
)

DOPC_m2 = Lipid(
    "DOPC",
    "NC3",
    "C1B",
    "C4B",
    12
)

