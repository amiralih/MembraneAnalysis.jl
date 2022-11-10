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

