"""
    pts, coeffs = initSource(mate::Material, mesh::Mesh, src::Source)

Source configuration. Returns source points `pts` and corresponding body force coefficients `coeffs`.
"""
function initSource(mate::Material, mesh::Mesh, src::Source)
    @unpack cp, c = mate
    @unpack spar, Ω = src
    @unpack type, dir, cmp = spar
    @unpack Nx, dx, X = mesh
    @bp
    
    pts = []
    coeffs = []
    if type == "Plane"
        # δ(x)
        # 1/Δx
        DIR = ([1,2,3] .== dir)
        nDIR = .!(DIR)
        pos = trunc(Int, 0.5*Nx[dir]) + 1
        imin = pos * DIR + nDIR
        imax = pos * DIR + (Nx .+ 1).*nDIR
        ind = [ [i,j,k] for i=imin[1]:imax[1], j=imin[2]:imax[2], k=imin[3]:imax[3] ]
        pts = [ ind[i] for i in eachindex(ind) ]
        force = zeros(Float64, 3)
        force[cmp] = 1.
        coeffs = [ force / dx[dir] for i in eachindex(ind) ]
    elseif type == "Cylinder" # polarised along z only (cmp = 3)
        ρ = 0.2
        pos = trunc.(Int, 0.5*Nx) .+ 1
        Xs = X[pos[1],pos[2],pos[3]]
        imin = pos .- trunc.(Int, ρ ./ dx) .- 1
        Imin = pos .- trunc.(Int, ρ/sqrt(2.) ./ dx)
        imax = pos .+ trunc.(Int, ρ ./ dx) .+ 1
        Imax = pos .+ trunc.(Int, ρ/sqrt(2.) ./ dx)
        imin = max.(imin, [1,1,1])
        imax = min.(imax, Nx .+ 1)
        Imin = max.(Imin, [1,1,1])
        Imax = min.(Imax, Nx .+ 1)
        ind = [ [i,j,k] for i=imin[1]:imax[1], j=imin[2]:imax[2], k=1:Nx[3]+1 ]
        Ind = [ [i,j,k] for i=Imin[1]:Imax[1], j=Imin[2]:Imax[2], k=1:Nx[3]+1 ]
        poi = [ ind[i] for i in eachindex(ind) ]
        Poi = [ Ind[i] for i in eachindex(Ind) ]
        pts = symdiff(poi, Poi)
        dists = [ (X[pt[1],pt[2],pt[3]] .- Xs) for pt in pts ]
        ang = [ atan(dists[pt][2], dists[pt][1]) for pt in eachindex(pts) ]
        diffs = [ abs.([dists[pt][1] - ρ*cos(ang[pt]), dists[pt][2] - ρ*sin(ang[pt])]) for pt in eachindex(pts) ]
        force = zeros(Float64, 3)
        force[cmp] = 1. / (dx[1]*dx[2])
        coeffs = [ (force * (diffs[pt][1] < 0.5*dx[1]) * (diffs[pt][2] < 0.5*dx[2])) for pt in eachindex(pts) ]
    elseif type == "Point"
        # coeffs = 1/√(2*π*σ^2)^3 exp(- X^T X / (2 σ^2)) with truncation and renormalisation
        R = 3.25
        σ = 0.018
        pos = trunc.(Int, 0.5*Nx) .+ 1
        Xs = X[pos[1],pos[2],pos[3]]
        imin = pos .- trunc.(Int, R*σ ./ dx) .- 1 # restriction |x| < 4 σ
        imax = pos .+ trunc.(Int, R*σ ./ dx) .+ 1
        imin = max.(imin, [1,1,1])
        imax = min.(imax, Nx .+ 1)
        ind = [ [i,j,k] for i=imin[1]:imax[1], j=imin[2]:imax[2], k=imin[3]:imax[3] ]
        pts = [ ind[i] for i in eachindex(ind) ]
        force = zeros(Float64, 3)
        cnorm = erf(R/sqrt(2)) - sqrt(2/pi) * R*exp(-0.5*R^2)
        force[cmp] = 1. / sqrt(2*π*σ^2)^3 / cnorm
        coeffs = [ (exp(-0.5*norm(X[pt[1],pt[2],pt[3]] .- Xs)^2/σ^2) * (norm(X[pt[1],pt[2],pt[3]] .- Xs)/σ < R) * force) for pt in pts ]
    else
        print("Unknown Source\n")
    end
    print("("*string(length(pts))*" source points set)\n")
    return pts, coeffs
end