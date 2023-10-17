"""
    xth, uth, res = compSol(mate::Material, mesh::Mesh, src::Source)

Computes analytical solutions via FFT.
"""
function compSol(mate::Material, mesh::Mesh, src::Source)
    @unpack type, dir, cmp = src.spar
    @unpack Ω, signal = src
    @unpack Nx, X, Qp, Nf, t = mesh
    @unpack ρ, μ, g, ω, c, Nv, ϵ = mate
    @bp
    
    if type == "Plane"
        # Slice
        DIR = ([1,2,3] .== dir)
        nDIR = .!(DIR)
        dirs = findall(==(0),DIR)
        poss = zeros(Int64,3)
        poss[dirs[1]] = trunc(Int, 0.5*Nx[dirs[1]]) + 1
        poss[dirs[2]] = trunc(Int, 0.5*Nx[dirs[2]]) + 1
        imin = poss .* nDIR + DIR
        imax = poss .* nDIR + (Nx .+ 1).*DIR
        xplot = [ X[i,j,k][dir] for i=imin[1]:imax[1], j=imin[2]:imax[2], k=imin[3]:imax[3] ]
        uplot = [ 0. for i=imin[1]:imax[1], j=imin[2]:imax[2], k=imin[3]:imax[3], f=1:Nf ]
        xth = copy( reshape(xplot, Nx[dir]+1) )
        uth = copy( reshape(uplot, (Nx[dir]+1, Nf)) )
        
        pos = trunc(Int, 0.5*Nx[dir]) + 1
        xsr = xth[pos*DIR + 1*nDIR]' * DIR
        force = zeros(3)
        force[cmp] = 1.
        dI = (cmp == dir)
        cI = c * sqrt(1. + (1/3. + 1/ϵ)*dI)
    
        if Nv <= 0
            for i=1:Nx[dir]+1
                uth[i,9+cmp] = 0.5/cI * signal(t - abs(xth[i]-xsr)/cI, Ω)
            end
        else
            # FFT analysis and synthesis
            # sampling
            T = 2*π/Ω
            Ne = floor(Int, 5e2)
            fe = Ne/T
            te = [0:Ne-1;]/fe
            se = [signal(te[i], Ω) for i=1:Ne]
            # zero-padding
            se = [se; zeros(nextpow(2,10*Ne)-Ne)]
            Nez = length(se)
            # fft
            she = rfft(se)/fe
            wez = 2*π*fe * rfftfreq(Nez)
            cI = vec(c * transpose(sqrt.( (1. .- sum(g .* ω ./ (ω .+ im*wez'), dims=1))*(1. + dI/3.) .+ dI/ϵ )))
            p = plan_irfft(she*fe, Nez; flags=FFTW.ESTIMATE, timelimit=Inf)
            fac = 0.5 * she*fe ./ cI
            shei = zeros(Nx[dir]+1)
            shhe = zeros(Nx[dir]+1)
            for i=1:Nx[dir]+1
                # ifft
                shei = fac .* exp.(im*wez .* (t .- abs(xth[i]-xsr) ./ cI))
                shhe = p * shei
                uth[i,9+cmp] += real( shhe[1] )
            end
        end
        
        return xth, uth, true
    else
        return [], [], false
        print("Unknown Source\n")
    end
end