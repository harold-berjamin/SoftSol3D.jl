"""
# ]
# (SoftSol3D) pkg> activate "docs"
# (docs) pkg> instantiate
#
# julia> include("docs/make.jl")
"""

# Packages
using Documenter, Dates
# ENV["JULIA_DEBUG"] = Documenter

# Modules
push!(LOAD_PATH,"../src/") # if path not loaded
using SoftSol3D

makedocs(;
    sitename = "SoftSol3D.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        edit_link = "https://github.com/harold-berjamin/SoftSol3D.jl.git",
        repolink = "https://github.com/harold-berjamin/SoftSol3D.jl.git"
    ),
    remotes = nothing,
    repo = Remotes.GitHub("harold-berjamin", "SoftSol3D.jl"),
    authors = "Harold Berjamin",
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "Manual" => "man.md",
    ]
)

deploydocs(;
    repo = "github.com/harold-berjamin/SoftSol3D.jl.git",
    target = "build",
    branch = "gh-pages",
    devbranch = "master",
    devurl = "dev",
    versions = ["stable" => "v^", "v#.#"],
)