"""
 Defines the data types using Parameters.jl.
"""
module ConfigPar

using Parameters, JuliaInterpreter, LinearAlgebra, SpecialFunctions
export Material, Mesh, Source
export SourcePar, initSource

# --------------------------------
# Physical domain

xmin = [-0.5, -0.1, -0.1] # [-0.5, -0.1, -0.1]
xmax = [ 0.5,  0.1,  0.1] # [ 0.5,  0.1,  0.1]
tmax = 0.20346395843291 # 0.20346395843291

# --------------------------------
# Constitutive behaviour

ρ = 1.000e3 # 1.000e3 (kg/m³)
μ = 2.684e3 # 2.684e3 (Pa)
β = 4.4 # 4.4

g = [0.0434, 0.0466, 0.2213] # [0.0434, 0.0466, 0.2213]
ω = 2*π * [10^1, 10^2, 10^3] # 2*π * [10^1, 10^2, 10^3]

# --------------------------------
# Computational domain

Nx = [100, 2, 2] # [100, 2, 2]
Nt = 0 # 0

# --------------------------------
# Computational method

flux = "LLF" # "LLF"
Co = 0.95 # 0.95 (Courant number)
ϵ = 0.9 # 0.9 (AC parameter)

# --------------------------------
# Loadings

type = "Plane" # "Plane", "SLine", "SPoint"
dir = 1 # orientation
cmp = 2 # component
Ω = 46.321609160460206 # 46.321609160460206 (angular frequency)

Sin2(t,Ω) = (sin(Ω*t) - 0.5*sin(2*Ω*t)) .* ((t>=0) - (t>=2*π/Ω))

# --------------------------------
# Data structures

"""
 Material parameters.
"""
@with_kw mutable struct Material
    ρ::Float64 = ρ; @assert ρ > 0.
    μ::Float64 = μ; @assert μ > 0.
    c::Float64 = sqrt(abs(μ/ρ))
    β::Float64 = β
    g::Array{Float64} = g
    ω::Array{Float64} = ω
    Nv::Int64 = length(g); @assert length(ω) == Nv
    ϵ::Float64 = ϵ; @assert ϵ > 0.
    cp::Float64 = sqrt(abs(μ/ρ*(4/3. + 1/ϵ))); @assert cp > 0.
end

"""
 Computational parameters.
"""
@with_kw mutable struct Mesh
    xmin::Array{Float64} = xmin; @assert length(xmin) == 3
    xmax::Array{Float64} = xmax; @assert length(xmax) == 3
    tmax::Float64 = tmax; @assert tmax >= 0.
    Nx::Array{Int64} = Nx; @assert length(Nx) == 3
    X::Array{Array{Float64,1},3} = [ xmin .+ ([i,j,k] .- 1) .* (xmax-xmin) ./ Nx for i=1:Nx[1]+1, j=1:Nx[2]+1, k=1:Nx[3]+1 ]
    Nt::Int64 = Nt; @assert Nt >= 0
    dx::Array{Float64} = (xmax-xmin) ./ Nx
    Co::Float64 = Co; @assert Co > 0.
    dt::Float64 = 0.1
    flux::String = flux
    t::Float64 = 0.
    n::Int64 = 0
    Nf::Int64 = 6*minimum(length, [g, ω]) + 12; @assert Nf > 0
    Qp::Array{Array{Float64,1},3} = [ [1.;0.;0.;0.;1.;0.;0.;0.;1.;zeros(Float64,Nf-9)] for i=1:Nx[1]+1, j=1:Nx[2]+1, k=1:Nx[3]+1 ]
    Fp05::Array{Array{Float64,1},3} = [ zeros(Nf) for i=1:Nx[1]+1, j=1:Nx[2]+1, k=1:Nx[3]+1 ]
end

"""
Source sub-parameters.
"""
@with_kw mutable struct SourcePar
    type::String = type
    dir::Int64 = dir; @assert (dir > 0) & (dir < 4)
    cmp::Int64 = cmp; @assert (cmp > 0) & (cmp < 4)
end

"""
 Source parameters (loading).
"""
@with_kw mutable struct Source
    spar::SourcePar = SourcePar()
    Ω::Float64 = Ω; @assert Ω > 0.
    signal::Function = Sin2
end

include("Sources.jl")

end # module