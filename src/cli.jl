# Command-line interface.

const CLI_USAGE = """
Usage:
  membrane-analysis [options] PARAMETERS.toml STRUCTURE.pdb TRAJECTORY [TRAJECTORY ...]

Analyse MD trajectories of a lipid membrane. PARAMETERS.toml defines the lipids, proteins and
analyses (see the documentation or examples/parameters.toml), STRUCTURE.pdb is the structure of
the system and the trajectories (any format readable by Chemfiles, e.g. .xtc, .dcd, .trr) are
given in chronological order.

Options:
  -o, --output DIR      output directory (overrides [output] directory in the parameters file)
  -n, --nprocs N        number of worker processes (default: CPUs allocated by SLURM/PBS/SGE,
                        otherwise all CPU threads)
      --resume          reuse per-trajectory results of an interrupted run in the output directory
      --chunk-frames N  split trajectories into work units of at most N frames
      --check           validate the inputs, print the analysis plan and exit
      --list-analyses   print the available analysis types and exit
  -h, --help            show this help and exit
      --version         show the package version and exit

Examples:
  membrane-analysis params.toml step5_charmm2gmx.pdb step7_{1..15}.xtc -o results -n 16
  julia --project=/path/to/MembraneAnalysis.jl /path/to/MembraneAnalysis.jl/bin/membrane-analysis.jl \\
        params.toml system.pdb traj_*.xtc
"""

"""
    parse_cli(args) -> NamedTuple

Parse command-line arguments. Throws [`ParameterError`](@ref) for invalid arguments.
"""
function parse_cli(args::AbstractVector{<:AbstractString})
    positional = String[]
    opts = Dict{Symbol,Any}(:output => nothing, :nprocs => nothing, :resume => false,
                            :check => false, :chunk_frames => nothing, :action => :run)
    i = 1
    function value_of(name)
        i < length(args) || throw(ParameterError("option $name needs a value"))
        i += 1
        return String(args[i])
    end
    function int_value(name, text)
        v = tryparse(Int, text)
        (v === nothing || v < 1) && throw(ParameterError("option $name needs a positive integer (got \"$text\")"))
        return v
    end
    while i <= length(args)
        a = String(args[i])
        if a == "-h" || a == "--help"
            opts[:action] = :help
        elseif a == "--version"
            opts[:action] = :version
        elseif a == "--list-analyses"
            opts[:action] = :list
        elseif a == "-o" || a == "--output"
            opts[:output] = value_of(a)
        elseif startswith(a, "--output=")
            opts[:output] = a[(length("--output=") + 1):end]
        elseif a == "-n" || a == "--nprocs"
            opts[:nprocs] = int_value(a, value_of(a))
        elseif startswith(a, "--nprocs=")
            opts[:nprocs] = int_value("--nprocs", a[(length("--nprocs=") + 1):end])
        elseif a == "--chunk-frames"
            opts[:chunk_frames] = int_value(a, value_of(a))
        elseif startswith(a, "--chunk-frames=")
            opts[:chunk_frames] = int_value("--chunk-frames", a[(length("--chunk-frames=") + 1):end])
        elseif a == "--resume"
            opts[:resume] = true
        elseif a == "--check"
            opts[:check] = true
        elseif a == "--"
            append!(positional, String.(args[(i + 1):end]))
            break
        elseif startswith(a, "-") && length(a) > 1
            throw(ParameterError("unknown option \"$a\" (see --help)"))
        else
            push!(positional, a)
        end
        i += 1
    end
    if opts[:action] === :run && length(positional) < 3
        throw(ParameterError("expected a parameters file, a structure file and at least one " *
                             "trajectory (see --help)"))
    end
    return (action=opts[:action], parameters=get(positional, 1, ""),
            structure=get(positional, 2, ""),
            trajectories=length(positional) >= 3 ? positional[3:end] : String[],
            output=opts[:output], nprocs=opts[:nprocs], resume=opts[:resume],
            check=opts[:check], chunk_frames=opts[:chunk_frames])
end

"""
    cli_main(args) -> Int

Entry point of the `membrane-analysis` command; returns the process exit code.
"""
function cli_main(args::AbstractVector{<:AbstractString})
    try
        o = parse_cli(args)
        if o.action === :help
            print(stdout, CLI_USAGE)
        elseif o.action === :version
            println(stdout, "MembraneAnalysis.jl v", PACKAGE_VERSION)
        elseif o.action === :list
            println(stdout, "Available analysis types:")
            for t in analysis_types()
                println(stdout, "  ", t)
            end
        else
            run_analysis(o.parameters, o.structure, o.trajectories; output_dir=o.output,
                         nprocs=o.nprocs, resume=o.resume, check=o.check,
                         chunk_frames=o.chunk_frames)
        end
        return 0
    catch err
        err isa InterruptException && (println(stderr, "interrupted"); return 130)
        if err isa ParameterError || err isa ArgumentError || err isa AnalysisError
            println(stderr, "error: ", err.msg)
        else
            println(stderr, "error: ", sprint(showerror, err))
            Base.show_backtrace(stderr, catch_backtrace())
            println(stderr)
        end
        return 1
    end
end

@static if isdefined(Base, Symbol("@main"))
    function (@main)(args)
        return cli_main(args)
    end
else
    """
        main(args) -> Int

    Entry point of the `membrane-analysis` command (same as [`cli_main`](@ref)).
    """
    main(args) = cli_main(args)
end
