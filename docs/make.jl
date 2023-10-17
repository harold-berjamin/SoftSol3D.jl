push!(LOAD_PATH,"../src/")

using Documenter, SoftSol3D

makedocs(
    sitename = "SoftSol3D.jl",
    modules  = [SoftSol3D],
    pages=[
        "Home" => "index.md"
    ]
)

deploydocs(;
    repo="github.com/harold-berjamin/SoftSol3D.jl",
)