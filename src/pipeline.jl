# Running a complete analysis: trajectories are split into work units that are processed in
# parallel (one frame in memory per process), per-trajectory time series are written to HDF5,
# and the time series of all trajectories are reduced to means and errors.

"""
    Job(system, analyses, surface, save_spectra)

Everything a worker process needs to process frames: sent once to every worker.
"""
struct Job
    system::MembraneSystem
    analyses::Vector{AbstractAnalysis}
    surface::SurfaceSettings
    save_spectra::Bool
end

struct TrajectoryInfo
    path::String
    n_frames::Int
end

struct WorkUnit
    trajectory::Int          # index of the trajectory in the input list
    path::String
    frames::UnitRange{Int}
end

struct UnitResult
    unit::WorkUnit
    sums::Vector{Matrix{Float64}}
    counts::Vector{Matrix{Float64}}
    box::Matrix{Float64}
    hq::Array{ComplexF32,3}
    tq::Array{ComplexF32,3}
    seconds::Float64
end

"""
    TrajectoryResult

Per-frame sums and counts of every analysis (`n_components × n_frames` matrices) and the box
lengths (`3 × n_frames`) of one trajectory.
"""
struct TrajectoryResult
    info::TrajectoryInfo
    sums::Vector{Matrix{Float64}}
    counts::Vector{Matrix{Float64}}
    box::Matrix{Float64}
end

"""
    AnalysisError(message)

An analysis run could not be completed (e.g. some trajectories failed); the message explains
what happened and how to continue.
"""
struct AnalysisError <: Exception
    msg::String
end
Base.showerror(io::IO, e::AnalysisError) = print(io, "AnalysisError: ", e.msg)

# ------------------------------------------------------------------------- frame loop

"""
    process_unit(job, unit) -> UnitResult

Single pass over the frames of one work unit: read a frame, centre it, compute the surface
(if needed) and let every analysis accumulate its per-frame values.
"""
function process_unit(job::Job, unit::WorkUnit)
    system = job.system
    n = length(unit.frames)
    sums = [zeros(n_components(a), n) for a in job.analyses]
    counts = [zeros(n_components(a), n) for a in job.analyses]
    box = zeros(3, n)
    need_surface = job.save_spectra || any(needs_surface, job.analyses)
    nx, ny = job.surface.grid_size
    ws = SurfaceWorkspace(job.surface, system)
    ws.leaflet .= system.static_leaflet
    hq = job.save_spectra ? zeros(ComplexF32, nx, ny, n) : zeros(ComplexF32, 0, 0, 0)
    tq = similar(hq)
    coords = Matrix{Float64}(undef, 3, system.n_atoms)
    t0 = time()
    with_trajectory(unit.path) do reader
        for (k, step) in enumerate(unit.frames)
            L = read_frame!(coords, reader, step)
            midplane = prepare_frame!(coords, L, system.tail_atoms;
                                      threshold=job.surface.center_threshold)
            need_surface && compute_surface!(ws, coords, L, midplane, system, job.surface)
            ctx = FrameContext(coords, L[1], L[2], L[3], midplane, ws.leaflet, ws.hq, ws.tq)
            for (a, analysis) in enumerate(job.analyses)
                accumulate_frame!(view(sums[a], :, k), view(counts[a], :, k), analysis, system, ctx)
            end
            box[1, k], box[2, k], box[3, k] = L
            if job.save_spectra
                hq[:, :, k] .= ws.hq
                tq[:, :, k] .= ws.tq
            end
        end
    end
    return UnitResult(unit, sums, counts, box, hq, tq, time() - t0)
end

"""
    foreach_frame(f, system, trajectory; surface=nothing, frames=nothing, center_threshold=20.0)

Memory-safe loop over the frames of a trajectory for custom analyses: calls `f(ctx, i)` with a
[`FrameContext`](@ref) for every frame `i` in `frames` (all frames by default). Coordinates are
wrapped and centred. With `surface::SurfaceSettings` the spectra and the dynamic lipid leaflets
are computed; otherwise `ctx.hq`/`ctx.tq` are empty and `ctx.leaflet` holds the leaflets of the
PDB structure. The `FrameContext` and its arrays are reused between frames — copy anything you
want to keep.
"""
function foreach_frame(f, system::MembraneSystem, path::AbstractString;
                       surface::Union{Nothing,SurfaceSettings}=nothing, frames=nothing,
                       center_threshold::Real=surface === nothing ? 20.0 : surface.center_threshold)
    coords = Matrix{Float64}(undef, 3, system.n_atoms)
    ws = surface === nothing ? nothing : SurfaceWorkspace(surface, system)
    leaflet = ws === nothing ? copy(system.static_leaflet) : ws.leaflet
    empty_spectrum = zeros(ComplexF64, 0, 0)
    with_trajectory(path) do reader
        for i in something(frames, 1:reader.n_frames)
            L = read_frame!(coords, reader, i)
            midplane = prepare_frame!(coords, L, system.tail_atoms; threshold=center_threshold)
            ws === nothing || compute_surface!(ws, coords, L, midplane, system, surface)
            ctx = FrameContext(coords, L[1], L[2], L[3], midplane, leaflet,
                               ws === nothing ? empty_spectrum : ws.hq,
                               ws === nothing ? empty_spectrum : ws.tq)
            f(ctx, i)
        end
    end
    return nothing
end

# --------------------------------------------------------------------- worker processes

const _WORKER_JOB = Ref{Any}(nothing)

_set_job!(job::Job) = (_WORKER_JOB[] = job; nothing)
_run_unit(unit::WorkUnit) = process_unit(_WORKER_JOB[]::Job, unit)

function _configure_worker()
    FFTW.set_num_threads(1)
    LinearAlgebra.BLAS.set_num_threads(1)
    return nothing
end

"""
    default_nprocs() -> Int

Number of CPUs available to this job: taken from the batch system (`SLURM_CPUS_PER_TASK`,
`SLURM_CPUS_ON_NODE`, `PBS_NP`, `NSLOTS`) if set, otherwise all CPU threads of the machine.
"""
function default_nprocs()
    for var in ("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE", "PBS_NP", "NSLOTS")
        v = tryparse(Int, get(ENV, var, ""))
        v !== nothing && v > 0 && return v
    end
    return Sys.CPU_THREADS
end

function _start_workers(n::Int)
    existing = filter(!=(1), workers())
    if !isempty(existing)
        _log("using $(length(existing)) existing worker processes")
        ws = existing
        added = Int[]
    else
        project = Base.active_project()
        flags = project === nothing ? `--threads=1` : `--threads=1 --project=$(dirname(project))`
        _log("starting $n worker processes")
        added = addprocs(n; exeflags=flags)
        ws = added
    end
    try
        Distributed.remotecall_eval(Main, ws, :(using MembraneAnalysis))
        Distributed.remotecall_eval(Main, ws, :(MembraneAnalysis._configure_worker()))
    catch
        isempty(added) || rmprocs(added)
        rethrow()
    end
    return ws, added
end

"""
    run_units(on_result, job, units, nprocs) -> failures

Process `units` with `nprocs` processes (serially in this process for `nprocs ≤ 1`), calling
`on_result(result)` in this process as soon as each unit is finished. Returns a list of
`(unit, error)` for failed units.
"""
function run_units(on_result, job::Job, units::Vector{WorkUnit}, nprocs::Int)
    failures = Tuple{WorkUnit,Any}[]
    no_workers = workers() == [1]
    if nprocs <= 1 || (length(units) <= 1 && no_workers)
        for u in units
            try
                on_result(process_unit(job, u))
            catch err
                err isa InterruptException && rethrow()
                push!(failures, (u, err))
                _log("FAILED: $(basename(u.path)) frames $(first(u.frames))–$(last(u.frames)): " *
                     sprint(showerror, err))
            end
        end
        return failures
    end
    ws, added = _start_workers(min(nprocs, length(units)))
    queue = collect(1:length(units))
    try
        @sync for w in ws
            @async remotecall_wait(_set_job!, w, job)
        end
        running = Ref(0)   # units currently being processed (they may be requeued)
        @sync for w in ws
            @async while true
                if isempty(queue)
                    running[] == 0 && break
                    sleep(0.2)
                    continue
                end
                i = popfirst!(queue)
                u = units[i]
                running[] += 1
                worker_died = false
                try
                    on_result(remotecall_fetch(_run_unit, w, u))
                catch err
                    if err isa ProcessExitedException || !(w in workers())
                        _log("worker $w died while processing $(basename(u.path)); the work unit is requeued")
                        push!(queue, i)
                        worker_died = true
                    else
                        push!(failures, (u, err))
                        _log("FAILED: $(basename(u.path)) frames $(first(u.frames))–$(last(u.frames)): " *
                             sprint(showerror, err))
                    end
                finally
                    running[] -= 1
                end
                worker_died && break
            end
        end
    catch
        empty!(queue)        # stop remaining tasks (e.g. after an interrupt)
        rethrow()
    finally
        isempty(added) || rmprocs(added)
    end
    for i in queue
        push!(failures, (units[i], ErrorException("no worker processes left")))
    end
    return failures
end

# ------------------------------------------------------------------------ work units

function make_units(infos::Vector{Tuple{Int,TrajectoryInfo}}, n_workers::Int;
                    chunk_frames::Union{Nothing,Int}=nothing)
    units = WorkUnit[]
    T = length(infos)
    for (k, info) in infos
        n = info.n_frames
        n_chunks = chunk_frames !== nothing ? cld(n, chunk_frames) :
                   T >= n_workers ? 1 : min(n, cld(n_workers, T))
        bounds = round.(Int, range(0, n; length=n_chunks + 1))
        for c in 1:n_chunks
            r = (bounds[c] + 1):bounds[c + 1]
            isempty(r) || push!(units, WorkUnit(k, info.path, r))
        end
    end
    return units
end

# ------------------------------------------------------------------ per-trajectory files

_timeseries_name(k::Int, path::AbstractString) = @sprintf("%03d_%s.h5", k, splitext(basename(path))[1])

function _fingerprint(cfg::AnalysisConfig, job::Job, info::TrajectoryInfo, structure::AbstractString)
    d = config_dict(cfg)
    delete!(d, "statistics")
    delete!(d, "output")
    delete!(d, "title")
    d["save_spectra"] = job.save_spectra
    d["grid_size"] = collect(job.surface.grid_size)
    d["structure"] = abspath(structure)
    d["trajectory"] = info.path
    d["n_frames"] = info.n_frames
    d["version"] = "$(PACKAGE_VERSION.major).$(PACKAGE_VERSION.minor)"
    d["layout"] = [vcat(analysis_id(a), component_labels(a)) for a in job.analyses]
    io = IOBuffer()
    TOML.print(io, d; sorted=true)
    return bytes2hex(SHA.sha1(take!(io)))
end

function write_trajectory_file(path::AbstractString, tr::TrajectoryResult, job::Job,
                               hq::Array{ComplexF32,3}, tq::Array{ComplexF32,3},
                               fingerprint::AbstractString)
    tmp = path * ".part"
    HDF5.h5open(tmp, "w") do f
        HDF5.write_attribute(f, "package_version", string(PACKAGE_VERSION))
        HDF5.write_attribute(f, "trajectory", tr.info.path)
        HDF5.write_attribute(f, "structure", job.system.pdb_file)
        HDF5.write_attribute(f, "n_frames", tr.info.n_frames)
        HDF5.write_attribute(f, "fingerprint", String(fingerprint))
        f["box"] = tr.box
        g = HDF5.create_group(f, "analyses")
        for (a, analysis) in enumerate(job.analyses)
            ga = HDF5.create_group(g, analysis_id(analysis))
            ga["values"] = tr.sums[a] ./ tr.counts[a]
            ga["counts"] = tr.counts[a]
            ga["sums"] = tr.sums[a]
            HDF5.write_attribute(ga, "labels", component_labels(analysis))
        end
        if !isempty(hq)
            gs = HDF5.create_group(f, "spectra")
            gs["hq"] = hq
            gs["tq"] = tq
            HDF5.write_attribute(gs, "grid_size", collect(job.surface.grid_size))
        end
    end
    mv(tmp, path; force=true)
    return path
end

function read_trajectory_file(path::AbstractString, job::Job, info::TrajectoryInfo,
                              fingerprint::AbstractString)
    isfile(path) || return nothing
    try
        return HDF5.h5open(path, "r") do f
            HDF5.read_attribute(f, "fingerprint") == fingerprint || return nothing
            HDF5.read_attribute(f, "n_frames") == info.n_frames || return nothing
            box = read(f["box"])
            sums = Matrix{Float64}[]
            counts = Matrix{Float64}[]
            for analysis in job.analyses
                ga = f["analyses"][analysis_id(analysis)]
                c = read(ga["counts"])
                if haskey(ga, "sums")
                    push!(sums, read(ga["sums"]))
                else
                    v = read(ga["values"])
                    push!(sums, ifelse.(c .== 0, 0.0, v .* c))
                end
                push!(counts, c)
            end
            TrajectoryResult(info, sums, counts, box)
        end
    catch err
        err isa InterruptException && rethrow()
        @warn "could not reuse $path; the trajectory will be processed again" exception = err
        return nothing
    end
end

# ---------------------------------------------------------------------------- main run

_log(msg::AbstractString) = (println(stdout, "[", Dates.format(Dates.now(), "HH:MM:SS"), "] ", msg); flush(stdout))

function _natural_key(path::AbstractString)
    return [(n = tryparse(Int, m.match); n === nothing ? (1, 0, String(m.match)) : (0, n, ""))
            for m in eachmatch(r"\d+|\D+", basename(path))]
end

function _scan_trajectories(paths::Vector{String}, system::MembraneSystem)
    max_frames = clamp(cld(500, length(paths)), 5, 50)
    infos = TrajectoryInfo[]
    sx = sy = sz = 0.0
    ns = 0
    for p in paths
        n, x, y, z, k = sample_box(p, system.n_atoms; max_frames=max_frames)
        push!(infos, TrajectoryInfo(p, n))
        sx += x; sy += y; sz += z; ns += k
    end
    return infos, (sx / ns, sy / ns, sz / ns)
end

function _describe_plan(cfg, system, infos, job, nprocs, outdir)
    lines = String[]
    push!(lines, "MembraneAnalysis.jl v$(PACKAGE_VERSION)" * (isempty(cfg.title) ? "" : " — " * cfg.title))
    push!(lines, sprint(show, MIME("text/plain"), system))
    push!(lines, "trajectories: $(length(infos)) files, $(sum(i.n_frames for i in infos)) frames")
    nx, ny = job.surface.grid_size
    s = job.surface
    surf = "surface: $(nx)×$(ny) grid, empty cells filled by $(s.fill)"
    s.mask_radius > 0 && (surf *= ", cells within $(s.mask_radius) Å of $(length(s.mask_atoms)) protein backbone atoms reconstructed")
    push!(lines, surf)
    push!(lines, "analyses:")
    for a in job.analyses
        push!(lines, "  - " * analysis_id(a) * " ($(n_components(a)) components)")
    end
    push!(lines, "statistics: blocking = $(cfg.blocking)" *
                 (cfg.blocking === :equal ? ", $(cfg.n_blocks) blocks" : ""))
    push!(lines, "processes: $nprocs; output directory: $outdir")
    return lines
end

"""
    run_analysis(parameters, structure, trajectories; output_dir=nothing, nprocs=nothing,
                 resume=false, check=false, chunk_frames=nothing) -> Dict

Run all analyses requested in the parameters file `parameters` (path or
[`AnalysisConfig`](@ref)) on the trajectories (in chronological order) of the system in the PDB
file `structure`.

Trajectories are processed in parallel by `nprocs` worker processes (default:
[`default_nprocs`](@ref)); long trajectories are split into chunks when there are fewer
trajectories than processes. Per-frame results of every trajectory are written to
`<output_dir>/timeseries/`, final means ± errors to `<output_dir>/*.dat` and scalar results to
`<output_dir>/summary.toml`, which is also returned as a dictionary.

- `output_dir`: overrides `[output] directory` of the parameters file.
- `resume = true`: reuse per-trajectory results from an earlier (interrupted) run with the same
  settings instead of processing those trajectories again.
- `check = true`: only validate the inputs and print the plan.
- `chunk_frames`: maximum number of frames per work unit (advanced).
"""
function run_analysis(parameters::AbstractString, structure::AbstractString,
                      trajectories::AbstractVector{<:AbstractString}; kwargs...)
    return run_analysis(load_parameters(parameters), structure, trajectories; kwargs...)
end

run_analysis(parameters::Union{AbstractString,AnalysisConfig}, structure::AbstractString,
             trajectory::AbstractString; kwargs...) =
    run_analysis(parameters, structure, [trajectory]; kwargs...)

function run_analysis(cfg::AnalysisConfig, structure::AbstractString,
                      trajectories::AbstractVector{<:AbstractString};
                      output_dir::Union{Nothing,AbstractString}=nothing,
                      nprocs::Union{Nothing,Integer}=nothing, resume::Bool=false,
                      check::Bool=false, chunk_frames::Union{Nothing,Integer}=nothing)
    isempty(trajectories) && throw(ParameterError("no trajectory files given"))
    isfile(structure) || throw(ParameterError("structure file not found: $structure"))
    paths = [abspath(String(t)) for t in trajectories]
    for p in paths
        isfile(p) || throw(ParameterError("trajectory file not found: $p"))
    end
    allunique(paths) || throw(ParameterError("the same trajectory file is given more than once"))
    if length(paths) > 1 && sort(paths; by=_natural_key) != paths
        @warn "The trajectories are not in natural (numerical) order. Block averaging assumes " *
              "that the given order is chronological."
    end
    outdir = abspath(something(output_dir, cfg.output_directory))
    existing_workers = filter(!=(1), workers())
    nprocs = Int(something(nprocs, isempty(existing_workers) ? default_nprocs() : length(existing_workers)))
    nprocs >= 1 || throw(ParameterError("nprocs must be at least 1"))
    chunk_frames === nothing || chunk_frames >= 1 ||
        throw(ParameterError("chunk_frames must be at least 1"))

    _log("reading structure $(structure)")
    system = MembraneSystem(structure, cfg.lipids, cfg.proteins)
    _log("scanning $(length(paths)) trajectories")
    infos, box_estimate = _scan_trajectories(paths, system)
    grid = something(cfg.grid_size, grid_size_for(box_estimate[1], box_estimate[2], cfg.grid_spacing))
    surface = surface_settings(cfg, system, grid)
    analyses = build_analyses(cfg, system; box_estimate=box_estimate, grid_size=grid)
    job = Job(system, analyses, surface, cfg.save_spectra)
    for line in _describe_plan(cfg, system, infos, job, nprocs, outdir)
        println(stdout, line)
    end
    check && return Dict{String,Any}()

    tsdir = joinpath(outdir, "timeseries")
    mkpath(tsdir)
    open(joinpath(outdir, "parameters_used.toml"), "w") do io
        TOML.print(io, config_dict(cfg); sorted=true)
    end

    T = length(infos)
    files = [joinpath(tsdir, _timeseries_name(k, info.path)) for (k, info) in enumerate(infos)]
    fingerprints = [_fingerprint(cfg, job, info, structure) for info in infos]
    results = Vector{Union{Nothing,TrajectoryResult}}(nothing, T)
    if resume
        for k in 1:T
            results[k] = read_trajectory_file(files[k], job, infos[k], fingerprints[k])
            results[k] === nothing || _log("reusing $(basename(files[k]))")
        end
    end
    todo = [(k, infos[k]) for k in 1:T if results[k] === nothing]
    n_parallel = isempty(existing_workers) ? nprocs : length(existing_workers)
    units = make_units(todo, n_parallel;
                       chunk_frames=chunk_frames === nothing ? nothing : Int(chunk_frames))
    expected = Dict{Int,Int}()
    for u in units
        expected[u.trajectory] = get(expected, u.trajectory, 0) + 1
    end
    pending = Dict{Int,Vector{UnitResult}}()
    n_done = Ref(0)

    function on_result(r::UnitResult)
        n_done[] += 1
        u = r.unit
        _log(@sprintf("[%d/%d] %s frames %d–%d done in %.1f s", n_done[], length(units),
                      basename(u.path), first(u.frames), last(u.frames), r.seconds))
        chunks = push!(get!(pending, u.trajectory, UnitResult[]), r)
        length(chunks) == expected[u.trajectory] || return nothing
        sort!(chunks; by=c -> first(c.unit.frames))
        k = u.trajectory
        na = length(job.analyses)
        tr = TrajectoryResult(infos[k],
                              [reduce(hcat, [c.sums[a] for c in chunks]) for a in 1:na],
                              [reduce(hcat, [c.counts[a] for c in chunks]) for a in 1:na],
                              reduce(hcat, [c.box for c in chunks]))
        hq = job.save_spectra ? cat((c.hq for c in chunks)...; dims=3) : zeros(ComplexF32, 0, 0, 0)
        tq = job.save_spectra ? cat((c.tq for c in chunks)...; dims=3) : zeros(ComplexF32, 0, 0, 0)
        write_trajectory_file(files[k], tr, job, hq, tq, fingerprints[k])
        results[k] = tr
        delete!(pending, k)
        return nothing
    end

    if !isempty(units)
        _log("processing $(length(units)) work units of $(length(todo)) trajectories")
        failures = run_units(on_result, job, units, nprocs)
        if !isempty(failures)
            failed = unique(basename(u.path) for (u, _) in failures)
            throw(AnalysisError(
                "$(length(failures)) work unit(s) failed (trajectories: $(join(failed, ", "))). " *
                "Completed trajectories are saved in $tsdir; fix the problem and rerun with " *
                "--resume (resume = true) to process only the remaining ones."))
        end
    end

    _log("computing statistics")
    summary = reduce_results(cfg, job, Vector{TrajectoryResult}(results), outdir, structure)
    _log("results written to $outdir")
    return summary
end

"""
    reduce_results(config, job, results, outdir, structure) -> Dict

Concatenate the per-frame results of all trajectories (in order), compute block statistics for
every analysis, write the result files and `summary.toml`.
"""
function reduce_results(cfg::AnalysisConfig, job::Job, results::Vector{TrajectoryResult},
                        outdir::AbstractString, structure::AbstractString)
    lengths = [r.info.n_frames for r in results]
    header = ["MembraneAnalysis.jl v$(PACKAGE_VERSION); structure: $(basename(structure)); " *
              "$(length(results)) trajectories, $(sum(lengths)) frames; blocking: $(cfg.blocking)"]
    out = OutputContext(outdir, header)
    summary = Dict{String,Any}()
    for (a, analysis) in enumerate(job.analyses)
        S = reduce(hcat, [r.sums[a] for r in results])
        C = reduce(hcat, [r.counts[a] for r in results])
        blocks = choose_blocks(cfg.blocking, S, C, lengths; n_blocks=cfg.n_blocks)
        bs = BlockStats(component_labels(analysis), S, C, blocks)
        res = write_results(analysis, bs, job.system, out)
        if !isempty(res)
            res["n_blocks"] = n_blocks(bs)
            summary[analysis_id(analysis)] = res
        end
    end
    summary["run"] = Dict{String,Any}(
        "package_version" => string(PACKAGE_VERSION),
        "date" => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"),
        "parameters" => cfg.source,
        "structure" => abspath(structure),
        "trajectories" => [r.info.path for r in results],
        "frames_per_trajectory" => lengths,
        "blocking" => String(cfg.blocking),
        "grid_size" => collect(job.surface.grid_size))
    open(joinpath(outdir, "summary.toml"), "w") do io
        TOML.print(io, summary; sorted=true)
    end
    return summary
end
