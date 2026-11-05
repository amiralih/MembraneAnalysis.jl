# Discretised membrane surface: leaflet heights on a grid, filling of empty/masked cells,
# dynamic leaflet assignment and the height (hq) and thickness (tq) fluctuation spectra.

"""
    SurfaceSettings(; grid_size, fill=:neighbors, mask_atoms=Int[], mask_radius=0.0,
                    center_threshold=20.0)

Settings of the membrane-surface calculation.

- `grid_size`: `(nx, ny)` number of grid cells.
- `fill`: how cells without reference atoms (and masked cells) are filled.
  `:neighbors` — sum of the four neighbouring cells (behaviour of earlier versions);
  `:min_curvature` — values that minimise the bending energy ∑ q⁴|h_q|² of the surface.
- `mask_atoms`, `mask_radius`: grid cells whose centre is closer than `mask_radius` (in the
  xy plane) to any of `mask_atoms` are treated as empty and reconstructed (requires
  `fill = :min_curvature`). Typically the backbone atoms of an embedded protein.
- `center_threshold`: see [`prepare_frame!`](@ref).
"""
struct SurfaceSettings
    grid_size::Tuple{Int,Int}
    fill::Symbol
    mask_atoms::Vector{Int}
    mask_radius::Float64
    center_threshold::Float64

    function SurfaceSettings(grid_size, fill, mask_atoms, mask_radius, center_threshold)
        nx, ny = grid_size
        (nx >= 4 && ny >= 4) || throw(ArgumentError("the surface grid must have at least 4×4 cells, got $(nx)×$(ny)"))
        fill in (:neighbors, :min_curvature) ||
            throw(ArgumentError("fill must be :neighbors or :min_curvature (got :$fill)"))
        if mask_radius > 0 && fill !== :min_curvature
            throw(ArgumentError("masking cells around atoms requires fill = :min_curvature"))
        end
        return new((Int(nx), Int(ny)), fill, collect(Int, mask_atoms), Float64(mask_radius),
                   Float64(center_threshold))
    end
end

SurfaceSettings(; grid_size, fill::Symbol=:neighbors, mask_atoms=Int[], mask_radius::Real=0.0,
                center_threshold::Real=20.0) =
    SurfaceSettings(grid_size, fill, mask_atoms, mask_radius, center_threshold)

"""
    grid_size_for(Lx, Ly, spacing) -> (nx, ny)

Number of grid cells so that cells are at most `spacing` long: `ceil(L / spacing)`.
"""
grid_size_for(Lx::Real, Ly::Real, spacing::Real) = (ceil(Int, Lx / spacing), ceil(Int, Ly / spacing))

"""
    SurfaceWorkspace(settings, system)
    SurfaceWorkspace(settings, n_molecules)

Preallocated buffers for [`compute_surface!`](@ref). One workspace is created per work unit (or
per [`foreach_frame`](@ref) call) and reused for every frame. After `compute_surface!` the fields `hq`, `tq` (spectra) and `leaflet`
(leaflet of every lipid molecule) hold the results for the current frame.
"""
struct SurfaceWorkspace
    nx::Int
    ny::Int
    zsum1::Matrix{Float64}
    zsum2::Matrix{Float64}
    n1::Matrix{Int}
    n2::Matrix{Int}
    h1::Matrix{Float64}
    h2::Matrix{Float64}
    hm::Matrix{Float64}
    protein_mask::BitMatrix
    mask::BitMatrix
    hq::Matrix{ComplexF64}
    tq::Matrix{ComplexF64}
    phase::Matrix{ComplexF64}   # shifts the FFT phase to cell centres
    q4::Matrix{Float64}         # bending-energy weights (i² + j²)² of the modes
    kernel::Matrix{Float64}     # real-space kernel of the bending energy, real(ifft(q4))
    leaflet::Vector{Int8}       # leaflet of every lipid molecule in the current frame
end

SurfaceWorkspace(settings::SurfaceSettings, system::MembraneSystem) =
    SurfaceWorkspace(settings, length(system.mol_head))

function SurfaceWorkspace(settings::SurfaceSettings, n_molecules::Integer)
    nx, ny = settings.grid_size
    q4 = bending_weights(nx, ny)
    kernel = real.(ifft(complex.(q4)))
    phase = [cis(-π * fft_index(i, nx) / nx - π * fft_index(j, ny) / ny) for i in 1:nx, j in 1:ny]
    return SurfaceWorkspace(nx, ny, zeros(nx, ny), zeros(nx, ny), zeros(Int, nx, ny),
                            zeros(Int, nx, ny), zeros(nx, ny), zeros(nx, ny), zeros(nx, ny),
                            falses(nx, ny), falses(nx, ny), zeros(ComplexF64, nx, ny),
                            zeros(ComplexF64, nx, ny), phase, q4, kernel,
                            zeros(Int8, n_molecules))
end

"""
    fft_index(i, n) -> Int

Signed frequency of the `i`-th (1-based) element of a length-`n` FFT (FFTW ordering:
0, 1, …, then negative frequencies).
"""
fft_index(i::Integer, n::Integer) = (k = i - 1; k >= cld(n, 2) ? k - n : k)

"""
    bending_weights(nx, ny) -> Matrix

Weights `(i² + j²)²` of the Fourier modes with `|i| ≤ nx÷2 - 1` and `|j| ≤ ny÷2 - 1` (zero for
all other modes), arranged in FFT order. ∑ weights .* |fft(h)|² is proportional to the bending
energy of the surface `h`.
"""
function bending_weights(nx::Int, ny::Int)
    w = zeros(nx, ny)
    mx = nx ÷ 2 - 1
    my = ny ÷ 2 - 1
    for i in -mx:mx, j in -my:my
        w[mod(i, nx) + 1, mod(j, ny) + 1] = (i^2 + j^2)^2
    end
    return w
end

# grid cell (1-based) of coordinate x in a box of length L divided into n cells
@inline _cell(x::Float64, n::Int, L::Float64) = clamp(floor(Int, x * n / L) + 1, 1, n)

function _bin_leaflet!(zsum, ns, coords, atoms, system, leaflet, which::Int8, Lx, Ly)
    fill!(zsum, 0.0)
    fill!(ns, 0)
    nx, ny = size(zsum)
    @inbounds for idx in atoms
        leaflet[system.atom_mol[idx]] == which || continue
        i = _cell(coords[1, idx], nx, Lx)
        j = _cell(coords[2, idx], ny, Ly)
        zsum[i, j] += coords[3, idx]
        ns[i, j] += 1
    end
    return nothing
end

"""
    fill_neighbors!(h, zsum, ns)

`h = zsum ./ ns` after giving every empty cell the summed height and count of its four
neighbours (cells are visited in order and filled in place, exactly as in earlier versions).
`zsum` and `ns` are modified.
"""
function fill_neighbors!(h::Matrix{Float64}, zsum::Matrix{Float64}, ns::Matrix{Int})
    nx, ny = size(zsum)
    for i in 1:nx, j in 1:ny
        if ns[i, j] == 0
            up    = (mod1(i - 1, nx), j)
            down  = (mod1(i + 1, nx), j)
            left  = (i, mod1(j - 1, ny))
            right = (i, mod1(j + 1, ny))
            zsum[i, j] = zsum[up...] + zsum[down...] + zsum[left...] + zsum[right...]
            ns[i, j] = ns[up...] + ns[down...] + ns[left...] + ns[right...]
        end
    end
    h .= zsum ./ ns
    return h
end

"""
    fill_min_curvature!(h, mask, ws)

Replace the values of `h` in the cells where `mask` is `true` by the values that minimise the
bending energy E = ∑_q q⁴ |h_q|² of the whole surface, keeping all other cells fixed.

Setting ∂E/∂h_k = 0 for the masked cells gives the linear system
`K[M, M] z = -(L h₀)[M]`, where `L = F⁻¹ q⁴ F` is the (circulant) bending operator, `K` its
real-space kernel, `h₀` the surface with zeros in the masked cells and `M` the masked cells.
This is algebraically identical to the normal equations built mode by mode in the
experimental `fluctuation_spectrum_hmin*` functions, but needs only one FFT pair per call.
"""
function fill_min_curvature!(h::Matrix{Float64}, mask::BitMatrix, ws::SurfaceWorkspace)
    cells = findall(mask)
    M = length(cells)
    M == 0 && return h
    nx, ny = size(h)
    for c in cells
        h[c] = 0.0
    end
    Lh = real.(ifft(ws.q4 .* fft(h)))
    A = Matrix{Float64}(undef, M, M)
    b = Vector{Float64}(undef, M)
    for k in 1:M
        ck = cells[k]
        b[k] = -Lh[ck]
        for l in 1:M
            cl = cells[l]
            A[k, l] = ws.kernel[mod(ck[1] - cl[1], nx) + 1, mod(ck[2] - cl[2], ny) + 1]
        end
    end
    z = _solve_symmetric(A, b)
    for k in 1:M
        h[cells[k]] = z[k]
    end
    return h
end

function _solve_symmetric(A::Matrix{Float64}, b::Vector{Float64})
    F = cholesky(Symmetric(A); check=false)
    issuccess(F) && return F \ b
    @warn "The surface reconstruction system is singular (too many masked cells?); " *
          "using the pseudo-inverse." maxlog = 1
    return pinv(A) * b
end

function _mark_protein_cells!(mask::BitMatrix, coords, settings::SurfaceSettings, Lx, Ly)
    fill!(mask, false)
    (settings.mask_radius > 0 && !isempty(settings.mask_atoms)) || return mask
    nx, ny = size(mask)
    r = settings.mask_radius
    @inbounds for i in 1:nx, j in 1:ny
        cx = (i - 0.5) * Lx / nx
        cy = (j - 0.5) * Ly / ny
        for idx in settings.mask_atoms
            dx = abs(cx - coords[1, idx]); dx = min(dx, Lx - dx)
            dy = abs(cy - coords[2, idx]); dy = min(dy, Ly - dy)
            if sqrt(dx^2 + dy^2) < r
                mask[i, j] = true
                break
            end
        end
    end
    return mask
end

"""
    compute_surface!(ws, coords, box, midplane, system, settings)

Calculate the membrane surface of one (prepared) frame:

1. bin the reference atoms of each leaflet (leaflets from the PDB structure) on the grid and
   fill empty/masked cells to obtain the midplane surface;
2. assign every lipid to the leaflet on whose side of the midplane surface its head atom is;
   if any lipid is assigned differently from the structure (e.g. after a flip-flop), steps 1–2
   are repeated once with the new leaflets;
3. bin the reference atoms again with these leaflets, fill empty/masked cells, subtract the
   midplane height and Fourier transform both leaflet surfaces.

Results: `ws.hq` (height spectrum, mean of the leaflets), `ws.tq` (thickness spectrum,
difference of the leaflets) and `ws.leaflet` (±1 for every lipid molecule). Spectra are
normalised by `Lx Ly / (nx ny)` and phase-shifted to the cell centres. The result depends only
on the frame itself, not on previously processed frames.
"""
function compute_surface!(ws::SurfaceWorkspace, coords::Matrix{Float64}, box::NTuple{3,Float64},
                          midplane::Float64, system::MembraneSystem, settings::SurfaceSettings)
    length(ws.leaflet) == length(system.mol_head) || throw(DimensionMismatch(
        "the workspace has $(length(ws.leaflet)) lipid molecules, the system $(length(system.mol_head))"))
    size(coords, 2) == system.n_atoms || throw(DimensionMismatch(
        "coords has $(size(coords, 2)) atoms, the system $(system.n_atoms)"))
    Lx, Ly, _ = box
    nx, ny = ws.nx, ws.ny
    atoms = system.surface_atoms
    minimise = settings.fill === :min_curvature
    if minimise
        _mark_protein_cells!(ws.protein_mask, coords, settings, Lx, Ly)
    end

    # leaflet assignment from the midplane surface (leaflets of the structure as a start)
    _assign_leaflets!(ws, coords, system, system.static_leaflet, Lx, Ly, minimise)
    if ws.leaflet != system.static_leaflet
        _assign_leaflets!(ws, coords, system, ws.leaflet, Lx, Ly, minimise)
    end

    # leaflet surfaces relative to the midplane
    _bin_leaflet!(ws.zsum1, ws.n1, coords, atoms, system, ws.leaflet, Int8(1), Lx, Ly)
    _bin_leaflet!(ws.zsum2, ws.n2, coords, atoms, system, ws.leaflet, Int8(-1), Lx, Ly)
    if minimise
        @. ws.h1 = ws.zsum1 / ws.n1 - midplane
        @. ws.h2 = ws.zsum2 / ws.n2 - midplane
        @. ws.mask = (ws.n1 == 0) | ws.protein_mask
        fill_min_curvature!(ws.h1, ws.mask, ws)
        @. ws.mask = (ws.n2 == 0) | ws.protein_mask
        fill_min_curvature!(ws.h2, ws.mask, ws)
    else
        fill_neighbors!(ws.h1, ws.zsum1, ws.n1)
        fill_neighbors!(ws.h2, ws.zsum2, ws.n2)
        ws.h1 .-= midplane
        ws.h2 .-= midplane
    end

    # spectra
    norm = Lx * Ly / (nx * ny)
    hq1 = fft(ws.h1) .* norm
    hq2 = fft(ws.h2) .* norm
    @. ws.hq = (hq1 + hq2) / 2 * ws.phase
    @. ws.tq = (hq1 - hq2) * ws.phase
    return ws
end

# Midplane surface (absolute height, in ws.hm) from the reference atoms with leaflets `source`,
# then assignment of every lipid to a leaflet (ws.leaflet). `source` may be ws.leaflet itself:
# it is only read before ws.leaflet is overwritten.
function _assign_leaflets!(ws::SurfaceWorkspace, coords::Matrix{Float64}, system::MembraneSystem,
                           source::Vector{Int8}, Lx::Float64, Ly::Float64, minimise::Bool)
    atoms = system.surface_atoms
    _bin_leaflet!(ws.zsum1, ws.n1, coords, atoms, system, source, Int8(1), Lx, Ly)
    _bin_leaflet!(ws.zsum2, ws.n2, coords, atoms, system, source, Int8(-1), Lx, Ly)
    if minimise
        @. ws.hm = (ws.zsum1 / ws.n1 + ws.zsum2 / ws.n2) / 2
        @. ws.mask = (ws.n1 == 0) | (ws.n2 == 0) | ws.protein_mask
        fill_min_curvature!(ws.hm, ws.mask, ws)
    else
        fill_neighbors!(ws.h1, ws.zsum1, ws.n1)
        fill_neighbors!(ws.h2, ws.zsum2, ws.n2)
        @. ws.hm = (ws.h1 + ws.h2) / 2
    end
    nx, ny = ws.nx, ws.ny
    @inbounds for m in eachindex(system.mol_head)
        a = system.mol_head[m]
        i = _cell(coords[1, a], nx, Lx)
        j = _cell(coords[2, a], ny, Ly)
        ws.leaflet[m] = coords[3, a] > ws.hm[i, j] ? Int8(1) : Int8(-1)
    end
    return nothing
end
