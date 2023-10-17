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
    elseif type == "Point"
        # coeffs = 1/√(2*π*σ^2)^3 exp(- X^T X / (2 σ^2)) with truncation and renormalisation
        R = 3.25
        σ = 0.018
        pos = trunc.(Int, 0.5*Nx) .+ 1
        Xs = X[pos]
        imin = pos .- trunc.(Int, R*σ ./ dx) .- 1 # restriction |x| < 4 σ
        imax = pos .+ trunc.(Int, R*σ ./ dx) .+ 1
        imin = max.(imin, [1,1,1])
        imax = min.(imax, Nx .+ 1)
        ind = [ [i,j,k] for i=imin[1]:imax[1], j=imin[2]:imax[2], k=imin[3]:imax[3] ]
        pts = [ ind[i] for i in eachindex(ind) ]
        force = zeros(Float64, 3)
        cnorm = erf(R/sqrt(2)) - sqrt(2/pi) * R*exp(-0.5*R^2)
        force[cmp] = 1. / sqrt(2*π*σ^2)^3 / cnorm
        coeffs = [ (exp(-0.5*norm(X[pt]-Xs)^2/σ^2) * (norm(X[pt]-Xs)/σ < R) * force) for pt in pts ]
    else
        print("Unknown Source\n")
    end
    print("("*string(length(pts))*" source points set)\n")
    return pts, coeffs
end