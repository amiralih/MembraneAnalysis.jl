# Analyses and output files

All quantities are averaged over frames (and over molecules/copies where applicable). The mean
of a quantity is ``\sum_\text{frames} S_f / \sum_\text{frames} n_f``, where ``S_f`` is the
sum of the values and ``n_f`` the number of values in frame ``f``. Standard errors (`sem`)
come from block averaging; derived quantities (``K_A``, ``k_c``, neutral surface, ``\Delta c``)
have block-jackknife errors (see [Methods](methods.md#Statistics)). Lengths are in Å, energies
in ``k_\mathrm{B}T``.

Signs: heights and curvatures of lipids are signed with respect to the lipid's leaflet
(positive towards the leaflet's surface; curvature positive when the surface bends towards the
leaflet). Proteins use the sign of their `leaflet` setting. Thickness values are not signed.

## `box`

Per frame: ``L_x, L_y, L_z, A = L_x L_y, A^2``.

| file | content |
|:-----|:--------|
| `box_dims.dat` | mean box lengths |
| `A.dat` | mean membrane area ``\langle A \rangle`` |
| `KA.dat` | area expansion modulus ``K_A = \langle A\rangle / (\langle A^2\rangle - \langle A\rangle^2)`` in ``k_\mathrm{B}T/\text{Å}^2`` |

## `undulation`

Per frame: ``\lvert h_{\mathbf q}\rvert^2`` and ``\lvert t_{\mathbf q}\rvert^2`` averaged over
the modes ``(n, m)``, ``0 \le n, m \le N``, with equal ``\sqrt{n^2 + m^2}``, where
``N = \max(\lfloor q_\text{max} L_x / 2\pi \rfloor, 4)``, limited to ``(\min(n_x, n_y) - 1) \div 2``
so that only resolved (non-aliased) modes are used.

| file | content |
|:-----|:--------|
| `hq2.dat` | ``q = 2\pi\sqrt{n^2+m^2}/\langle L_x\rangle`` and ``\langle \lvert h_{\mathbf q}\rvert^2 \rangle`` |
| `tq2.dat` | same for the thickness modes |
| `kc_hq2.dat` | bending modulus from a least-squares fit of ``q^4 \langle \lvert h_{\mathbf q}\rvert^2\rangle`` to the chosen model of [Terzi & Deserno, *J. Chem. Phys.* 147, 084702 (2017)](https://doi.org/10.1063/1.4990404): `HC` ``= p_1`` (eq. 36), `MNK` ``= p_1 + p_2 q^2`` (eq. 37), `TD` ``= (p_1 + p_2 q^2)/(1 - q^2/p_3)`` (eq. 47); ``k_c = \langle A\rangle / p_1`` |

The spectra assume a square box for the definition of ``q``.

## `lipid_heights`, `protein_heights`

Mean height above the bilayer midplane of every height atom of each lipid species
(`<LIPID>_zs.dat`: atom, mean, sem) or of the backbone atom(s) of every protein residue
(`<PROTEIN>_zs.dat`: residue, mean, sem).

## `lipid_sampled_curvature`, `protein_sampled_curvature`

Mean curvature per mode sampled at the position of each height atom (lipids, `<LIPID>_cs.dat`)
or residue backbone atom (proteins, `<PROTEIN>_cs.dat`), from the modes with
``q_\text{min} < |\mathbf q| < q_\text{max}``:

```math
c(\mathbf r) = \frac{1}{n/2}\,\mathrm{Re} \sum_{q_\text{min} < |\mathbf q| < q_\text{max}}
               \frac{q^2 h_{\mathbf q}}{L_x L_y}\, e^{i \mathbf q \cdot \mathbf r},
```

where ``n`` is the number of modes in the sum.

## `lipid_curvature_spectrum`, `protein_curvature_spectrum`

The same quantity in `n_bins` shells of ``|\mathbf q|`` between `q_min` and `q_max`, evaluated
at the reference atom of each lipid (`<LIPID>_cqs.dat`) or at the backbone atom(s) of the
protein residue `residue` (`<PROTEIN>_cqs.dat`). Columns: bin centre ``q``, mean, sem. Bins that
contain no Fourier mode in any frame give `NaN`. Because the box fluctuates, a mode close to a
bin edge can fall into one bin in some frames and into the neighbouring bin in others; each bin
is then averaged over the frames in which it contains modes (v0.3 reported `NaN` for both
bins). The number of samples per bin and frame is stored in `timeseries/*.h5` (`counts`).

## `lipid_thickness_spectrum`, `protein_thickness_spectrum`

As the curvature spectra, but from the thickness spectrum ``t_{\mathbf q}`` and without the
factor ``q^2`` (`<LIPID>_tqs.dat`, `<PROTEIN>_tqs.dat`).

## Directional analyses

With `direction = "x"` (`"y"`), only modes with wave vectors along x (y), i.e. ``(n, 0)``
(``(0, m)``), are used. Default file names get the suffix `_X` (`_Y`), e.g. `PEP_cqs_X.dat`.

## `tcb`

Transverse curvature bias analysis ([Sapp, Beaven & Sodt, *Phys. Rev. E* 103, 042413 (2021)](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.103.042413)).
For each lipid species, the mean sampled curvature ``\langle c\rangle`` of every height atom is
fitted linearly against the atom's mean height ``z`` (atoms with ``z <`` `z_cutoff`):
``\langle c \rangle = a_\ell + b_\ell z``. With area fractions ``w_\ell``,
``b = \sum_\ell w_\ell b_\ell`` and ``a = \sum_\ell w_\ell a_\ell``:

| file | content |
|:-----|:--------|
| `TCB.dat` | bending modulus ``k_c = -2/(\langle A\rangle b)`` and neutral surface ``-a/b`` |
| `TCB_lipids.dat` | per lipid: weight, ``a_\ell``, ``b_\ell`` and the curvature bias ``\Delta c_\ell = a_\ell - a`` |
| `<LIPID>_TCB_profile.dat` | ``z`` and ``\langle c\rangle`` of every atom (for plotting the TCB curves) |

A TCB plot can be made from the profile files, for example with Plots.jl:

```julia
using DelimitedFiles, Plots
p = plot(xlabel="z (Å)", ylabel="⟨c⟩(z) (1/Å)")
for lipid in ("DOPC", "CHL1")
    d = readdlm("results/$(lipid)_TCB_profile.dat"; comments=true)
    scatter!(p, Float64.(d[:, 2]), Float64.(d[:, 4]); yerror=Float64.(d[:, 5]), label=lipid)
end
savefig(p, "TCB.pdf")
```
