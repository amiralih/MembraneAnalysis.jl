using Documenter
using MembraneAnalysis

makedocs(
    sitename = "MembraneAnalysis.jl",
    modules = [MembraneAnalysis],
    pages = [
        "Introduction" => "index.md",
        "Getting started" => "getting_started.md",
        "Parameters file" => "parameters.md",
        "Analyses and output files" => "analyses.md",
        "Methods" => "methods.md",
        "Julia API" => "api.md",
        "Migrating from v0.3" => "migration.md",
    ],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true",
                             assets = String[]),
    warnonly = [:missing_docs, :cross_references],
)

deploydocs(
    repo = "github.com/amiralih/MembraneAnalysis.jl.git",
    devbranch = "main",
)
