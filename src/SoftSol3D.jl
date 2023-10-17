module SoftSol3D

using Parameters, LinearAlgebra, JuliaInterpreter, FFTW, StaticArrays
export Run, Slice, compSol

include("ConfigPar.jl")
using .ConfigPar: Material, Mesh, Source
using .ConfigPar: SourcePar, initSource
export Material, Mesh, Source, LocalMaxSpeeds
export SourcePar, initSource

"""
    Run(mate::Material, mesh::Mesh, src::Source)

Executes the main loop for the parameters specified in the argument. The numerical solution 
at the final time `data.t` is stored in `data.Qp`. This data vector of size ``6N+13`` has the following structure:

``{\\bf Q} = \\left(F_{11}, F_{12}, F_{13}, \\dots , F_{33}, v_1, v_2, v_3, [S_1^\\text{v}]_{11}, \\dots , [S_N^\\text{v}]_{33}\\right)^\\top``

where ``F_{ij}`` are the deformation gradient components, ``v_i`` are the velocity components, and ``[S_\\ell^\\text{v}]_{ij}`` are the memory variables.
"""
function Run(mate::Material, mesh::Mesh, src::Source)
    
    # extract data
    @unpack ρ, μ, c, β, g, ω, Nv, ϵ, cp = mate
    @unpack tmax, Nx, X, Nt, dx, Co, dt, flux = mesh
    @unpack t, n, Nf, Qp, Fp05 = mesh
    @unpack Ω, signal = src
    @bp

    print(string("Δx: ", dx,"\n"))

    # source configuration
    pts, coeffs3 = initSource(mate, mesh, src)
    coeffs = [ [zeros(Float64,9); coeffs3[c]; zeros(Float64,Nf-12)] for c in eachindex(coeffs3) ]
    
    # time step initialisation
    cmax = MaxSpeeds(Qp, mate)
    dt = Co*minimum(dx ./ cmax)
    @pack! mesh = dt

    # numerical flux configuration
    fluxNum, fluxPhy, locVelo = getFlux(mate, mesh)
    
    while (t+dt <= tmax) && ((n+1 < Nt) || (Nt <= 0))
        # relaxation step using matrix exp
        ExpRel(Qp, mesh, mate)
        
        # propagation step
        # dimensional splitting (spatial dimension dir)
        for dir=1:3
            # compute flux
            fluxNum(Qp, Fp05, dir, mate, mesh, fluxPhy, locVelo)
            # update fields
            FVDiff(Fp05, Qp, dir, mesh)
            # boundaries (interior points)
            GhostExtra(Qp, dir, mesh)
        end

        # relaxation step using matrix exp
        ExpRel(Qp, mesh, mate)

        # source
        for pt in eachindex(pts)
            ind = pts[pt]
            Qp[ind[1],ind[2],ind[3]] .+= dt * signal(t+dt,Ω) * coeffs[pt] 
        end

        # update and increment
        t += dt
        n += 1
        @pack! mesh = t, n

        # maximum sound speeds
        cmax = MaxSpeeds(Qp, mate)
        dt = Co*minimum(dx ./ cmax)
        @pack! mesh = dt
    end

    # output
    @pack! mesh = t, n, Qp
    return 0
end

"""
    LocalMaxSpeeds(Q::Array{Float64,1}, mate::Material, cmaxloc::Array{Float64,1})

Returns the local maximum characteristic speeds along x, y, z at the state vector Q with the
parameters of the argument.
"""
function LocalMaxSpeeds(Q::Array{Float64,1}, mate::Material, cmaxloc::Array{Float64,1})
    @unpack ρ, μ, ϵ, β, Nv, cp = mate
    @bp
    
    K = μ/ϵ
    for i=1:3
        cmaxloc[i] = cp
    end
    FT = SMatrix{3,3}(@view(Q[1:9]))
    F = FT'
    J = det(F)
    if J <= 0.
        throw(DomainError(J, "Jacobian must be positive"))
    end
    Jm23 = J^(-2/3.)
    I1 = F ⋅ F
    I1b = I1 * Jm23
    W1b = 0.5*μ*(1. + 2/3. * β*(I1b-3.))
    a0 = 2.0*W1b * Jm23
    a1 = 4/3. * Jm23^2 * μ*β
    a2 = -4/3. * Jm23
    q = K*J^2 + 5/9. * 2.0*I1b*W1b
    if Nv>0
        Sv = zeros(Float64,3,3)
        sv = @view(Q[13:18])
        for ℓ = 2:Nv
            sv += @view(Q[13+6*(ℓ-1):18+6*(ℓ-1)])
        end
        ind = [1 2 3; 2 4 5; 3 5 6]
        for it in eachindex(Sv)
            Sv[it] = sv[ind[it]]
        end
        τv = F * Sv * FT
        q -= 5/9. * Jm23 * tr(τv)
    end

    FmT = inv(FT)
    N = zeros(Float64,3)
    m = @SVector zeros(Float64,3)
    md = @SVector zeros(Float64,3)
    for dir=1:3
        N[dir] = 1.
        m = FmT * N
        md = F * N # m^*
        ma = 2.0*W1b*md
        md -= I1/3. * m
        if Nv>0
            ma -= τv*m
            a0 -= dot(m,τv,m) * Jm23
        end
        Mat = Symmetric( q * m*m' + a0 * I + a1 * md*md' + 0.5*a2*(m*ma' + ma*m') )
        Matval = eigvals(Mat)
        if any(x -> x >= 0., Matval)
            cmaxloc[dir] = max(cmaxloc[dir], sqrt(maximum(Matval)/ρ))
        else
            throw(DomainError(Matval, "Eigenvalues must be real-valued non-negative"))
        end
        N[dir] = 0.
    end
end

"""
    LocalMaxSpeedsLin(Q::Array{Float64,1}, mate::Material, cmaxloc::Array{Float64,1})

LocalMaxSpeeds clone for the linear case
"""
function LocalMaxSpeedsLin(Q::Array{Float64,1}, mate::Material, cmaxloc::Array{Float64,1})
    @unpack cp = mate
    @bp

    for i=1:3
        cmaxloc[i] = cp
    end
end

"""
    MaxSpeeds(Q::Array{Array{Float64,1},3}, mate::Material)::Tuple{Float64, Float64, Float64}

Returns the global maximum characteristic speeds along x, y, z.
"""
function MaxSpeeds(Q::Array{Array{Float64,1},3}, mate::Material)::Tuple{Float64, Float64, Float64}
    @unpack cp, β = mate
    @bp

    cmax1 = cp
    cmax2 = cp
    cmax3 = cp
    cmaxloc = zeros(Float64,3)
    if abs(β) > eps(Float64)
        Nxp1 = size(Q)
        for i=1:Nxp1[1], j=1:Nxp1[2], k=1:Nxp1[3]
            LocalMaxSpeeds(Q[i,j,k], mate, cmaxloc)
            cmax1 = max(cmax1, cmaxloc[1])
            cmax2 = max(cmax2, cmaxloc[2])
            cmax3 = max(cmax3, cmaxloc[3])
        end
    end
    return cmax1, cmax2, cmax3
end

"""
    PhysFlux(Q::Array{Float64,1}, dir::Int64, mate::Material, flx::Array{Float64,1})

Returns the physical flux vector along direction `dir`.
"""
function PhysFlux(Q::Array{Float64,1}, dir::Int64, mate::Material, flx::Array{Float64,1})
    @unpack ρ, μ, β, Nv, ϵ = mate
    @bp

    K = μ/ϵ
    FT = SMatrix{3,3}(@view(Q[1:9]))
    F = FT'
    J = det(F)
    Jm23 = J^(-2/3.)
    if J <= 0.
        throw(DomainError(J, "Jacobian must be positive"))
    end
    Fm = inv(F)
    FmT = Fm'
    Fb = J^(-1/3) * F
    I1 = Fb ⋅ Fb
    W1 = 0.5*μ*(1. + 2/3. * β*(I1-3.))
    S = (K*J*(J-1) - 2/3. *W1*I1) * Fm*FmT + Jm23 * 2.0*W1*I
    if Nv>0
        Sv = zeros(Float64,3,3)
        sv = @view(Q[13:18])
        for ℓ = 2:Nv
            sv += @view(Q[13+6*(ℓ-1):18+6*(ℓ-1)])
        end
        ind = [1 2 3; 2 4 5; 3 5 6]
        for it in eachindex(Sv)
            Sv[it] = sv[ind[it]]
        end
        S -= Jm23 * (Sv - 1/3. * tr(F*Sv*FT) * Fm*FmT)
    end
    P = F * S
    flx[dir]   = -Q[10]
    flx[3+dir] = -Q[11]
    flx[6+dir] = -Q[12]
    flx[10] = -P[1,dir]/ρ
    flx[11] = -P[2,dir]/ρ
    flx[12] = -P[3,dir]/ρ
    return 0
end

"""
    PhysFluxLin(Q::Array{Float64,1}, dir::Int64, mate::Material, flx::Array{Float64,1})

Returns the linearised physical flux vector along direction `dir`.
"""
function PhysFluxLin(Q::Array{Float64,1}, dir::Int64, mate::Material, flx::Array{Float64,1})
    @unpack ρ, μ, β, Nv, ϵ = mate
    @bp

    F = SMatrix{3,3}(@view(Q[1:9]))'
    P = μ*(F+F' - 2.0*I) + (1/ϵ-2/3.)*μ * (tr(F)-3.0) * I
    if Nv>0
        Sv = zeros(Float64,3,3)
        sv = @view(Q[13:18])
        for ℓ = 2:Nv
            sv += @view(Q[13+6*(ℓ-1):18+6*(ℓ-1)])
        end
        ind = [1 2 3; 2 4 5; 3 5 6]
        for it in eachindex(Sv)
            Sv[it] = sv[ind[it]]
        end
        P -= Sv
    end
    flx[dir]   = -Q[10]
    flx[3+dir] = -Q[11]
    flx[6+dir] = -Q[12]
    flx[10] = -P[1,dir]/ρ
    flx[11] = -P[2,dir]/ρ
    flx[12] = -P[3,dir]/ρ
    return 0
end

"""
    FVDiff(Fp05::Array{Array{Float64,1},3}, Qp::Array{Array{Float64,1},3}, dir::Int64, mesh::Mesh)

Updates `Q` by adding the difference of numerical fluxes `Fp05` along `dir`.
"""
function FVDiff(Fp05::Array{Array{Float64,1},3}, Qp::Array{Array{Float64,1},3}, dir::Int64,
    mesh::Mesh)
    @unpack Nx, Nf, dx, dt, flux = mesh
    @bp

    ind = ([1,2,3] .== dir)
    if flux == "LLF"
        imin = 1 .+ ind
        imax = Nx .+ 1 .- ind
    elseif flux == "MUSCL"
        imin = 1 .+ 2*ind
        imax = Nx .+ 1 .- 2*ind
    end

    for i=imin[1]:imax[1], j=imin[2]:imax[2], k=imin[3]:imax[3]
        fp05 = Fp05[i,j,k]
        fm05 = Fp05[i-ind[1],j-ind[2],k-ind[3]]
        Qp[i,j,k] -= dt/dx[dir] * (fp05 - fm05)
    end
    return 0
end

"""
    GhostExtra(Qp::Array{Array{Float64,1},3}, dir::Int64, mesh::Mesh)

Extrapolates to ghost cell values along direction `dir`.
"""
function GhostExtra(Qp::Array{Array{Float64,1},3}, dir::Int64, mesh::Mesh)
    @unpack Nx, flux = mesh
    @bp

    ind = ([1,2,3] .== dir)
    nind = .!(ind)
    if flux == "LLF"
        imin = ind + nind
        imax = ind + nind .* (Nx .+ 1)
        icop = 2
    elseif flux == "MUSCL"
        imin =   ind + nind
        imax = 2*ind + nind .* (Nx .+ 1)
        icop = 3
    end

    if (Nx[dir]+2-icop < icop)
        throw(DomainError(icop, "Boundaries must be thicker. Enlarge computational domain."))
    end

    for i=imin[1]:imax[1], j=imin[2]:imax[2], k=imin[3]:imax[3]
        il, jl, kl = icop*ind .+ nind .* [i,j,k]
        Qp[i,j,k] .= Qp[il, jl, kl]
        ir, jr, kr = (Nx[dir]+2-icop).*ind .+ nind .* [i,j,k]
        ic, jc, kc = (Nx .+ 2 .- [i,j,k]).*ind + nind .* [i,j,k]
        Qp[ic, jc, kc] .= Qp[ir, jr, kr]
    end

    return 0
end

"""
    ExpRel(Qp::Array{Array{Float64,1},3}, mesh::Mesh, mate::Material)

Performs the relaxation step using matrix exponentials.
"""
function ExpRel(Qp::Array{Array{Float64,1},3}, mesh::Mesh, mate::Material)
    @unpack dt, Nx = mesh
    @unpack μ, β, g, ω, Nv = mate
    @bp
    
    if Nv>0
        if abs(β) > eps(Float64)
            for i=1:Nx[1]+1, j=1:Nx[2]+1, k=1:Nx[3]+1
                F = SMatrix{3,3}(@view(Qp[i,j,k][1:9]))'
                J = det(F)
                Fb = J^(-1/3.) * F
                Fm = inv(F)
                I1 = Fb ⋅ Fb
                W1 = 0.5*μ*(1. + 2/3. * β*(I1-3.))
                SeD = -2/3. * I1*W1 * Fm*Fm' + 2.0*W1*I
                seD = [@view(SeD[1,1:3]); @view(SeD[2,2:3]); SeD[3,3]]
                for ℓ=1:Nv
                    expL = exp(-0.5*ω[ℓ]*dt)
                    Qp[i,j,k][13+6*(ℓ-1):18+6*(ℓ-1)] = expL * @view(Qp[i,j,k][13+6*(ℓ-1):18+6*(ℓ-1)]) + (1. - expL)*g[ℓ]*seD
                end
            end
        else
            for i=1:Nx[1]+1, j=1:Nx[2]+1, k=1:Nx[3]+1
                F = reshape(@view(Qp[i,j,k][1:9]),3,3)' # no use of SMatrix since no det, inv, etc.
                SeD = μ*(F+F' - 2.0*I - 2/3. * tr(F-I)*I)
                seD = [@view(SeD[1,1:3]); @view(SeD[2,2:3]); SeD[3,3]]
                for ℓ=1:Nv
                    expL = exp(-0.5*ω[ℓ]*dt)
                    Qp[i,j,k][13+6*(ℓ-1):18+6*(ℓ-1)] = expL * @view(Qp[i,j,k][13+6*(ℓ-1):18+6*(ℓ-1)]) + (1. - expL)*g[ℓ]*seD
                end
            end
        end
    end
    return 0
end

"""
    x, u = Slice(dir::Int64, mesh::Mesh, data::Data)

Returns data `data.Qp` sliced along direction `dir` in output `u`. Spatial coordinates `mesh.X` in output `x` undergo the same operation.
"""
function Slice(dir::Int64, mesh::Mesh)
    @unpack Nx, X, Qp, Nf = mesh
    @bp
    
    DIR = ([1,2,3] .== dir)
    nDIR = .!(DIR)
    dirs = findall(==(0),DIR)
    poss = zeros(Int64,3)
    poss[dirs[1]] = trunc(Int, 0.5*Nx[dirs[1]]) + 1
    poss[dirs[2]] = trunc(Int, 0.5*Nx[dirs[2]]) + 1
    imin = poss .* nDIR + DIR
    imax = poss .* nDIR + (Nx .+ 1).*DIR
    xplot = [ X[i,j,k][dir] for i=imin[1]:imax[1], j=imin[2]:imax[2], k=imin[3]:imax[3] ]
    uplot = [ Qp[i,j,k][f] for i=imin[1]:imax[1], j=imin[2]:imax[2], k=imin[3]:imax[3], f=1:Nf ]
    return reshape(xplot, Nx[dir]+1), reshape(uplot, (Nx[dir]+1, Nf))
end

include("Schemes.jl")
include("CompTheo.jl")

end # module