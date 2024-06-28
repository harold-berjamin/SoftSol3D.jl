The reference documentation.

### Getting started

Open Julia's REPL, navigate to the downloaded SoftSol3D folder located at `/home/JuliaUser/Downloads/SoftSol3D` (see documentation homepage for the installation). Then activate and instantiate the package. To do so, type `]` in the REPL and run:

```julia
(SoftSol3D) pkg> activate .
(SoftSol3D) pkg> instantiate
```

### Running a simulation

Configure the parameters in `src/main.jl`, located at `/home/JuliaUser/Downloads/SoftSol3D`. In the REPL, run:

```julia
julia> include("src/main.jl")
```

### Running examples

In the REPL, run:

```julia
julia> include("test/main_Cauchy1D.jl")
```