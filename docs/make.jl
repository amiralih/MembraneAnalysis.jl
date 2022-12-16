using Documenter
using MembraneAnalysis

push!(LOAD_PATH,"../src/")
makedocs(
    sitename="MembraneAnalysis.jl",
    pages = [
             "Index" => "index.md",
             "Curvature Analysis" => "curvature_analysis.md",
             "Another Page" => "anotherPage.md",
             "Tutorial" => "anotherPage.md",
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
