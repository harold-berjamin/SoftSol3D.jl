"""
# ]
# (SoftSol3D) pkg> activate .
# (SoftSol3D) pkg> instantiate
#
# julia> include("test/main_Cauchy1D.jl")
"""

# Load 'Dev' packages
using Revise, BenchmarkTools

# Load 'Run' packages
using CSV, DataFrames, Dates
using Plots, WriteVTK
gr() # Plots backend (gr, pyplot, plotly, inspectdr, unicodeplots)

# Load modules
using Parameters, JuliaInterpreter, SoftSol3D
print("Dependencies loaded\n")

# Main function (no declaration in global scope)
function main()

    # Configuration
    mate = Material(
        ρ = 1.0e3, # 1.0e3 (kg/m³)
        μ = 2.684e3, # 2.684e3 (Pa)
        β = 0., # 4.4
        g = [], # [0.0434, 0.0466, 0.2213]
        ω = 2*π * [], # 2*π * [10^1, 10^2, 10^3]
        ϵ = 0.9 # 0.9, 4*(1/50)^0.3 (AC parameter)
    )
    mesh = Mesh(
        xmin = [-0.5, -0.1, -0.1], # [-0.5, -0.1, -0.1]
        xmax = [ 0.5,  0.1,  0.1], # [ 0.5,  0.1,  0.1]
        tmax = 0.18, # 0.18
        Nx = [1000, 4, 4], # [100, 2, 2], [100, 4, 4]
        Nt = 0, # 0
        Co = 0.95, # 0.95 (Courant number)
        flux = "MUSCL", # "LLF", "MUSCL"
        Nf = 6*minimum(length, [mate.g, mate.ω]) + 12 # 6*minimum(length, [mate.g, mate.ω]) + 12
    )
    Sin2(t::Float64, Ω::Float64) = 1e0 * (sin(Ω*t) - 0.5*sin(2*Ω*t)) .* ((t>=0) - (t>=2*π/Ω))
    src = Source(
        spar = SourcePar(
            type = "None", # "None", "Plane", "Cylinder", "Point"
            dir = 1, # 1
            cmp = 2 # 2
        ),
        Ω = 46.3, # 46.3
        signal = Sin2 # Sin2
    )
    fpara = 9 + src.spar.cmp # field for paraview display
    ana = true # analytical computations
    print("Configuration defined\n")

    # Initialisation
    @unpack c = mate
    @unpack Nx, X, Qp = mesh
    @unpack Ω, signal = src
    @bp

    for i=1:Nx[1]+1, j=1:Nx[2]+1, k=1:Nx[3]+1
        Qp[i,j,k][4] = -0.5*signal(-X[i,1,1][1]/c,Ω)/c^2
        Qp[i,j,k][11] = 0.5*signal(-X[i,1,1][1]/c,Ω)/c
    end

    # Iterations
    print("Run start (" * Dates.format(now(), "HH:MM") * ")\n")
    @time Run(mate, mesh, src)
    print("Run end (" * Dates.format(now(), "HH:MM") * ")\n")

    # Plot centered slice
    xplot, uplot = Slice(src.spar.dir, mesh)
    hp1 = plot(xplot, uplot[:,10:12], title = string("t = ", mesh.t), markershape = :circle, markersize=3)
    display(hp1)
    xlabel!("space (m)")
    ylabel!("velocities (m/s)")

    # Write text output
    df = DataFrame([xplot uplot], :auto)
    CSV.write("output/Result.csv", df)
    print("Results saved in output/Result.csv\n")

    # Write VTK output
    if fpara > 0
        xpara = mesh.xmin[1]:mesh.dx[1]:mesh.xmax[1]
        ypara = mesh.xmin[2]:mesh.dx[2]:mesh.xmax[2]
        zpara = mesh.xmin[3]:mesh.dx[3]:mesh.xmax[3]
        Nxpara = [length(xpara), length(ypara), length(zpara)]
        upara = [ mesh.Qp[i,j,k][fpara] for i=1:Nxpara[1], j=1:Nxpara[2], k=1:Nxpara[3] ]
        vtk_grid("output/Results", xpara, ypara, zpara) do vtk
            vtk["output"] = upara
        end
        print("Results saved in output/Results.vti\n")
    end

    # Analytical solution
    if ana
        # xth, uth, res = compSol(mate, mesh, src)
        
        res = false
        jy = trunc(Int, 0.5*Nx[2]) + 1
        kz = trunc(Int, 0.5*Nx[3]) + 1
        xp = [ X[i,jy,kz][1] for i=1:Nx[1]+1 ]
        up = [ 0. for i=1:Nx[1]+1, f=1:mesh.Nf ]
        xth = copy( reshape(xp, Nx[1]+1) )
        uth = copy( reshape(up, (Nx[1]+1, mesh.Nf)) )
        if mate.Nv <= 0
            for i=1:Nx[1]+1
                uth[i,11] = 0.5/c * signal(mesh.t - xth[i]/c, Ω)
            end
            res = true
        end

        if res
            hp2 = plot!(hp1, xth, uth[:,10:12], color="black")
            display(hp2)

            # Write text output
            dfth = DataFrame([xth uth], :auto)
            CSV.write("output/ResultTh.csv", dfth)
            print("Results saved in output/ResultTh.csv\n")
        end
    end

    print(" ")
end

main()