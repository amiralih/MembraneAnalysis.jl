using Documenter
using MembraneAnalysis

push!(LOAD_PATH,"../src/")
makedocs(
    sitename="MembraneAnalysis.jl",
    pages = [
             "Introduction" => "index.md",
             "Fluctuation Analysis" => "fluctuation_analysis.md",
             "Curvature Analysis" => "curvature_analysis.md",
             "Thickness Analysis" => "thickness_analysis.md",
             "Lipids" => "lipids.md",
             "Density Analysis" => "density_analysis.md",
             "Utilities" => "utilities.md",
            ],
    format = Documenter.HTML(prettyurls = false)
)
# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/amiralih/MembraneAnalysis.jl.git",
    devbranch = "main"
)
