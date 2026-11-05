# Memory-safe trajectory access.
#
# Earlier versions called `Chemfiles.read_step(traj, i)` for every frame. Each call allocates a
# new `Chemfiles.Frame` whose memory lives in the C++ library and is released only when Julia's
# garbage collector finalises the small wrapper object. The collector does not see the C++
# memory, so it has no reason to run promptly, and the memory of many frames can accumulate.
# Here a single `Frame` is allocated per open trajectory and reused with `Chemfiles.read_step!`,
# and coordinates are copied into one preallocated Julia matrix.

"""
    TrajectoryReader(path)

Open a trajectory file (any format supported by Chemfiles) for reading. Frames are read one at
a time into a single reused buffer, so memory use does not grow with the number of frames.
Close the reader with `close(reader)` or use [`with_trajectory`](@ref).
"""
mutable struct TrajectoryReader
    path::String
    trajectory::Chemfiles.Trajectory
    frame::Chemfiles.Frame
    n_frames::Int
end

function TrajectoryReader(path::AbstractString)
    isfile(path) || throw(ArgumentError("trajectory file not found: $path"))
    trajectory = Chemfiles.Trajectory(String(path))
    n_frames = Int(Chemfiles.size(trajectory))
    return TrajectoryReader(String(path), trajectory, Chemfiles.Frame(), n_frames)
end

Base.close(r::TrajectoryReader) = close(r.trajectory)
Base.length(r::TrajectoryReader) = r.n_frames

"""
    with_trajectory(f, path)

Open `path`, call `f(reader)` and close the trajectory afterwards, even if `f` throws.
"""
function with_trajectory(f, path::AbstractString)
    reader = TrajectoryReader(path)
    try
        return f(reader)
    finally
        close(reader)
    end
end

"""
    read_frame!(coords, reader, i) -> (Lx, Ly, Lz)

Read frame `i` (1-based) of `reader` into the preallocated `3 × N` matrix `coords` and return
the box lengths. Only rectangular boxes are supported.
"""
function read_frame!(coords::Matrix{Float64}, r::TrajectoryReader, i::Integer)
    1 <= i <= r.n_frames || throw(BoundsError(r, i))
    Chemfiles.read_step!(r.trajectory, i - 1, r.frame)
    positions = Chemfiles.positions(r.frame)
    if size(positions) != size(coords)
        throw(ArgumentError("frame $i of $(r.path) has $(size(positions, 2)) atoms but the " *
                            "structure has $(size(coords, 2)) atoms"))
    end
    copyto!(coords, positions)
    cell = Chemfiles.UnitCell(r.frame)
    angles = Chemfiles.angles(cell)
    if any(a -> abs(a - 90) > 1e-3, angles)
        throw(ArgumentError("frame $i of $(r.path) has a non-rectangular box (angles $(angles)); " *
                            "only rectangular boxes are supported"))
    end
    L = Chemfiles.lengths(cell)
    all(>(0), L) || throw(ArgumentError("frame $i of $(r.path) has no periodic box (box lengths $(L))"))
    return (Float64(L[1]), Float64(L[2]), Float64(L[3]))
end

"""
    trajectory_length(path) -> Int

Number of frames in a trajectory file.
"""
trajectory_length(path::AbstractString) = with_trajectory(r -> r.n_frames, path)

"""
    sample_box(path, n_atoms; max_frames=50) -> (n_frames, Lx_sum, Ly_sum, Lz_sum, n_sampled)

Read up to `max_frames` evenly spaced frames of a trajectory and return the number of frames
in the file together with sums of the box lengths over the sampled frames.
"""
function sample_box(path::AbstractString, n_atoms::Integer; max_frames::Integer=50)
    return with_trajectory(path) do r
        n = r.n_frames
        n == 0 && throw(ArgumentError("trajectory $path contains no frames"))
        steps = n <= max_frames ? (1:n) :
                unique(round.(Int, range(1, n; length=max_frames)))
        coords = Matrix{Float64}(undef, 3, n_atoms)
        sums = [0.0, 0.0, 0.0]
        for i in steps
            L = read_frame!(coords, r, i)
            sums .+= L
        end
        (n, sums[1], sums[2], sums[3], length(steps))
    end
end

"""
    prepare_frame!(coords, box, tail_atoms; threshold=20.0) -> midplane

Wrap all atoms into the box and shift them along z so that the bilayer lies in the middle of
the box (same procedure as in earlier versions):

1. while the standard deviation of the tail-atom heights exceeds `threshold` (the bilayer is
   split across the periodic boundary), shift all atoms by 20 Å along z;
2. if the midplane (mean tail-atom height) is further than 5 % of the box height from the box
   centre, shift all atoms so that it is at the centre.

Returns the height of the bilayer midplane.
"""
function prepare_frame!(coords::Matrix{Float64}, box::NTuple{3,Float64},
                        tail_atoms::Vector{Int}; threshold::Real=20.0)
    Lx, Ly, Lz = box
    @inbounds for i in axes(coords, 2)
        coords[1, i] = mod(coords[1, i], Lx)
        coords[2, i] = mod(coords[2, i], Ly)
        coords[3, i] = mod(coords[3, i], Lz)
    end
    tails_z = view(coords, 3, tail_atoms)
    max_shifts = ceil(Int, Lz / 20.0) + 1
    n_shifts = 0
    while std(tails_z) > threshold
        n_shifts += 1
        n_shifts > max_shifts && error(
            "could not centre the bilayer: the standard deviation of the tail-atom heights stays " *
            "above $threshold Å for every shift along z. Is the membrane flat?")
        _shift_z!(coords, 20.0, Lz)
    end
    midplane = mean(tails_z)
    Δz = Lz / 2 - midplane
    if abs(Δz) > 0.05 * Lz
        _shift_z!(coords, Δz, Lz)
        midplane = mean(tails_z)
    end
    return midplane
end

function _shift_z!(coords::Matrix{Float64}, Δz::Real, Lz::Real)
    @inbounds for i in axes(coords, 2)
        coords[3, i] = mod(coords[3, i] + Δz, Lz)
    end
    return coords
end
