# Methods

## Frame preparation

Every frame is processed independently:

1. All atoms are wrapped into the (rectangular) box.
2. While the standard deviation of the heights of the lipid tail atoms exceeds
   `center_threshold` (20 Å), the bilayer is split by the periodic boundary in z and all atoms
   are shifted by 20 Å in z. An error is raised if no shift brings the standard deviation below
   the threshold (e.g. for strongly curved membranes).
3. The midplane ``z_0`` is the mean height of all tail atoms. If it is further than 5 % of
   ``L_z`` from the box centre, all atoms are shifted so that it is at the centre.

## Membrane surface

The xy plane is divided into ``n_x \times n_y`` cells. The grid is the same for all frames and
trajectories: ``n = \lceil \langle L \rangle / \ell \rceil`` with the grid spacing ``\ell`` and
the mean box length ``\langle L\rangle`` estimated from evenly spaced frames of every trajectory
(up to 50 per trajectory, about 500 in total; or set explicitly with `grid_size`). In every frame, the cell of an atom at ``x`` is
``\lfloor x\, n_x / L_x \rfloor + 1`` with the current box length.

1. **Midplane surface.** The reference atoms of each leaflet (leaflets as in the PDB
   structure) are binned and the mean height per cell is computed for both leaflets; the
   midplane surface is their average. Cells without atoms (and masked cells, see below) are
   filled as described below.
2. **Leaflets.** Every lipid is assigned to the upper leaflet if its head atom lies above the
   midplane surface at its position and to the lower leaflet otherwise. If any lipid ends up in
   a different leaflet than in the structure (flip-flop), steps 1–2 are repeated once with the
   new leaflets. These dynamic leaflets are used by all analyses. The result depends only on the
   frame itself, so it does not depend on how trajectories are split between processes.
3. **Leaflet surfaces.** The reference atoms are binned again with the dynamic leaflets, giving
   the leaflet heights ``h_1, h_2`` relative to ``z_0`` (empty/masked cells filled).
4. **Spectra.** ``h_{\mathbf q}^{(i)} = \frac{L_x L_y}{n_x n_y}\, \mathrm{FFT}(h_i)``, phase-shifted
   by ``e^{-i\pi(k_x/n_x + k_y/n_y)}`` so that the modes refer to cell centres; the height
   spectrum is ``h_{\mathbf q} = (h^{(1)}_{\mathbf q} + h^{(2)}_{\mathbf q})/2`` and the
   thickness spectrum ``t_{\mathbf q} = h^{(1)}_{\mathbf q} - h^{(2)}_{\mathbf q}``.

### Filling empty cells: `fill = "neighbors"`

Cells are visited in order; an empty cell receives the summed heights and counts of its four
neighbours (including cells filled earlier in the same pass). This is the method of v0.3.

### Filling empty and masked cells: `fill = "min_curvature"`

Cells without atoms and — with `mask_radius > 0` — cells whose centre is closer than
`mask_radius` to a backbone atom of a masked protein are treated as unknown. Their values are
chosen to minimise the bending energy of the surface,

```math
E[h] = \sum_{|i| \le n_x/2 - 1,\ |j| \le n_y/2 - 1} (i^2 + j^2)^2 \, |\hat h_{ij}|^2 ,
```

with all other cells fixed. Writing ``h = h_0 + \sum_{k \in M} z_k e_k`` (``h_0``: the surface
with zeros in the set ``M`` of unknown cells, ``e_k``: unit surface of cell ``k``) and setting
``\partial E / \partial z_k = 0`` gives the linear system

```math
\sum_{l \in M} K(\mathbf r_k - \mathbf r_l)\, z_l = -(\mathcal L h_0)(\mathbf r_k), \qquad k \in M,
```

where ``\mathcal L = \mathrm{FFT}^{-1} w\, \mathrm{FFT}`` is the circulant bending operator with
weights ``w_{ij} = (i^2+j^2)^2`` and ``K = \mathrm{FFT}^{-1}(w)`` its real-space kernel. The
system is symmetric positive (semi-)definite and is solved by a Cholesky factorisation
(pseudo-inverse if singular). This is algebraically identical to the normal equations built mode
by mode in the experimental `fluctuation_spectrum_hmin*` functions (the test suite checks this),
but requires a single FFT pair per reconstruction instead of one FFT per unknown cell.

In the first step the midplane surface is reconstructed; in step 3 each leaflet surface is
reconstructed separately (its own empty cells plus the masked cells).

## Fourier evaluation at points

Curvature and thickness at a point ``\mathbf r`` are evaluated from the modes in a shell
``q_\text{min} < |\mathbf q| < q_\text{max}`` (``|i| \le n_x/2 - 1``, ``|j| \le n_y/2 - 1``):

```math
c(\mathbf r) = \frac{2}{n}\,\mathrm{Re}\sum_{\mathbf q} \frac{q^2 h_{\mathbf q}}{L_x L_y}
               e^{i\mathbf q\cdot\mathbf r}, \qquad
t(\mathbf r) = \frac{2}{n}\,\mathrm{Re}\sum_{\mathbf q} \frac{t_{\mathbf q}}{L_x L_y}
               e^{i\mathbf q\cdot\mathbf r},
```

with ``n`` the number of modes in the shell (``n/2`` independent modes, since
``h_{-\mathbf q} = h_{\mathbf q}^*``), i.e. values *per mode*. Directional analyses use only the
modes ``(i, 0)`` or ``(0, j)``.

## Statistics

Every analysis accumulates, for each frame ``f`` and component ``c`` (a q bin, an atom, …), the
sum of values ``S_{cf}`` and their number ``n_{cf}``. The frames of all trajectories are
concatenated in the given order.

- **Mean:** ``\bar x_c = \sum_f S_{cf} / \sum_f n_{cf}``.
- **Blocks:** with `blocking = "equal"`, ``B`` (`n_blocks`) contiguous blocks of
  ``\lfloor N/B \rfloor`` frames (remaining frames are used for the mean but not for the error);
  with `"trajectory"`, one block per trajectory file; with `"auto"`, blocks of the length found
  by the automated blocking method of
  [Jonsson, *Phys. Rev. E* 98, 043304 (2018)](https://doi.org/10.1103/PhysRevE.98.043304)
  (Flyvbjerg–Petersen blocking with a χ² test), taking the largest length over the components
  of the analysis.
- **Standard error:** ``\mathrm{sem}_c = s(\bar x_{c,b}) / \sqrt{B}``, where ``\bar x_{c,b}`` are
  the block means and ``s`` their standard deviation.
- **Derived quantities** ``\theta = f(\bar x)`` (``K_A``, ``k_c``, …) get block-jackknife errors:
  ``\theta_b = f(\bar x^{(-b)})`` with block ``b`` left out and
  ``\sigma_\theta^2 = \frac{B-1}{B}\sum_b (\theta_b - \bar\theta)^2``.

Compared with v0.3, which formed blocks from the list of all individual values, blocks now
always consist of whole frames; means are unchanged, standard errors of quantities averaged over
many molecules differ somewhat.

## Memory and parallelism

Each process opens its trajectory once and reads frames with `Chemfiles.read_step!` into a
single reused `Chemfiles.Frame`; coordinates are copied into one preallocated matrix. (Earlier
versions created a new `Frame` for every step. Its memory is owned by the C++ library and is
released only when Julia's garbage collector finalises the small wrapper object; the collector
does not see the C++ memory, so the memory of many frames could accumulate.) Per-frame results are a few numbers per analysis and component.
