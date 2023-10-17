push!(LOAD_PATH,"../src/") # if path not loaded

# Packages
using Documenter, Dates

# Modules
using SoftSol3D

makedocs(
    sitename = "SoftSol3D.jl",
    authors = "Harold Berjamin",
    modules  = [SoftSol3D],
    pages=[
        "Home" => "index.md",
        "Manual" => "man.md",
        "API" => "api.md"
    ]
)

deploydocs(;
    repo="github.com/harold-berjamin/SoftSol3D.jl",
)