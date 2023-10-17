"""
    getFlux(mate::Material, mesh::Mesh)::Function

Returns the desired numerical flux procedure.
"""
function getFlux(mate::Material, mesh::Mesh)::Tuple{Function, Function, Function}
    @unpack β = mate
    @unpack flux = mesh
    @bp

    fluxPhy = PhysFlux
    locVelo = LocalMaxSpeeds
    fluxNum = FluxLLF

    if flux == "MUSCL"
        fluxNum = FluxMUSCL
    end
    if abs(β)<eps(Float64)
        fluxPhy = PhysFluxLin
        locVelo = LocalMaxSpeedsLin
    end
    return fluxNum, fluxPhy, locVelo
end

"""
    FluxLLF(Q::Array{Array{Float64,1},3}, Fp05::Array{Array{Float64,1},3}, dir::Int64, mate::Material, mesh::Mesh, fluxPhy, locVelo)

Inserts the Rusanov fluxes along direction `dir` in `Fp05`.
"""
function FluxLLF(Q::Array{Array{Float64,1},3}, Fp05::Array{Array{Float64,1},3}, dir::Int64,
    mate::Material, mesh::Mesh, fluxPhy, locVelo)
    @unpack Nx, Nf = mesh

    ind = ([1,2,3] .== dir)
    imin = [1,1,1]
    imax = Nx .+ 1 .- ind

    fp0 = zeros(Float64,Nf)
    fp1 = zeros(Float64,Nf)
    cp0 = zeros(Float64,3)
    cp1 = zeros(Float64,3)

    for i=imin[1]:imax[1], j=imin[2]:imax[2], k=imin[3]:imax[3]
        Qp0 = Q[i,j,k]
        Qp1 = Q[i+ind[1],j+ind[2],k+ind[3]]
        fluxPhy(Qp0,dir,mate,fp0)
        fluxPhy(Qp1,dir,mate,fp1)
        locVelo(Qp0,mate,cp0)
        locVelo(Qp1,mate,cp1)
        cploc = max(cp0[dir],cp1[dir])
        Fp05[i,j,k] = 0.5*(fp0 + fp1) - 0.5*cploc*(Qp1 - Qp0)
    end
    return 0
end

"""
    zSlope(a::Array{Float64,1}, b::Array{Float64,1})::Array{Float64,1}

Returns zero slopes (first-order scheme).
"""
function zSlope(a::Array{Float64,1}, b::Array{Float64,1})::Array{Float64,1}
    return 0 .* a
end

"""
    stSlope(a::Array{Float64,1}, b::Array{Float64,1})::Array{Float64,1}

Returns the slope corresponding to classical schemes (Beam-Warming, Lax-Wendroff, Fromm = default).
Selection via line comments in source code.
"""
function stSlope(a::Array{Float64,1}, b::Array{Float64,1})::Array{Float64,1}
    # return a # Beam-Warming
    # return b # Lax-Wendroff
    return 0.5*(a+b) # Fromm (default)
end

"""
    mmSlope(a::Array{Float64,1}, b::Array{Float64,1})::Array{Float64,1}

Returns the minmod slopes.
"""
function mmSlope(a::Array{Float64,1}, b::Array{Float64,1})::Array{Float64,1}
    return 0.5*(sign.(a)+sign.(b)) .* min.(abs.(a),abs.(b))
end

"""
    MCSlope(a::Array{Float64,1}, b::Array{Float64,1})::Array{Float64,1}

Returns the monotonised central-difference slopes.
"""
function MCSlope(a::Array{Float64,1}, b::Array{Float64,1})::Array{Float64,1}
    return 0.5*(sign.(a)+sign.(b)) .* min.(2*abs.(a),2*abs.(b),0.5*abs.(a+b))
end

"""
    FluxMUSCL(Q::Array{Array{Float64,1},3}, Fp05::Array{Array{Float64,1},3}, dir::Int64, mate::Material, mesh::Mesh, fluxPhy, locVelo)

Inserts the MUSCL fluxes along direction `dir` in `Fp05` for linear solids.
"""
function FluxMUSCL(Q::Array{Array{Float64,1},3}, Fp05::Array{Array{Float64,1},3}, dir::Int64,
    mate::Material, mesh::Mesh, fluxPhy, locVelo)
    @unpack Nx, Nf, dx, dt = mesh
    @bp

    ind = ([1,2,3] .== dir)
    imin = [1,1,1] .+ ind
    imax = Nx .+ 1 .- 2*ind

    fp0Rp = zeros(Float64,Nf)
    fp0Rm = zeros(Float64,Nf)
    fp1Lp = zeros(Float64,Nf)
    fp1Lm = zeros(Float64,Nf)
    fp0R = zeros(Float64,Nf)
    fp1L = zeros(Float64,Nf)
    cp0R = zeros(Float64,3)
    cp1L = zeros(Float64,3)

    for i=imin[1]:imax[1], j=imin[2]:imax[2], k=imin[3]:imax[3]
        Qp0 = Q[i,j,k]
        Qp1 = Q[i+ind[1],j+ind[2],k+ind[3]]
        Qp2 = Q[i+2*ind[1],j+2*ind[2],k+2*ind[3]]
        Qm1 = Q[i-ind[1],j-ind[2],k-ind[3]]
        Δp0 = mmSlope(Qp0-Qm1, Qp1-Qp0)
        Δp1 = mmSlope(Qp1-Qp0, Qp2-Qp1)
        fluxPhy(Qp0+0.5*Δp0,dir,mate,fp0Rp)
        fluxPhy(Qp0-0.5*Δp0,dir,mate,fp0Rm)
        Qp0R = Qp0 + 0.5*Δp0 + 0.5*dt/dx[dir]*(fp0Rm-fp0Rp)
        fluxPhy(Qp1+0.5*Δp1,dir,mate,fp1Lp)
        fluxPhy(Qp1-0.5*Δp1,dir,mate,fp1Lm)
        Qp1L = Qp1 - 0.5*Δp1 + 0.5*dt/dx[dir]*(fp1Lm-fp1Lp)
        fluxPhy(Qp0R,dir,mate,fp0R)
        fluxPhy(Qp1L,dir,mate,fp1L)
        locVelo(Qp0R,mate,cp0R)
        locVelo(Qp1L,mate,cp1L)
        cploc = max(cp0R[dir],cp1L[dir])
        Fp05[i,j,k] = 0.5*(fp0R + fp1L) - 0.5*cploc*(Qp1L - Qp0R)
    end
    return 0
end