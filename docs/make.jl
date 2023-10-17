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

makedocs(
    sitename = "SoftSol3D.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        edit_link = nothing,
        repolink = "github.com/harold-berjamin/SoftSol3D.jl.git"
    ),
    remotes = nothing,
    authors = "Harold Berjamin",
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "Manual" => "man.md",
    ],
    doctest = true
)

deploydocs(;
    repo = "github.com/harold-berjamin/SoftSol3D.jl.git",
    devbranch = "master",
    push_preview = true
)