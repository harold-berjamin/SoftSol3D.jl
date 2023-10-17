"""
# ]
# (SoftSol3D) pkg> activate .
# (SoftSol3D) pkg> instantiate
#
# julia> include("test/main_Shear3D.jl")
"""

# Load 'Dev' packages
using Revise, BenchmarkTools

# Load 'Run' packages
using CSV, DataFrames, Dates
using Plots, WriteVTK
gr() # Plots backend (gr, pyplot, plotly, inspectdr, unicodeplots)

# Load modules
using SoftSol3D
print("Dependencies loaded\n")

# Main function (no declaration in global scope)
function main()

    # Configuration
    mate = Material(
        ρ = 1.0e3, # 1.0e3 (kg/m³)
        μ = 2.684e3, # 2.684e3 (Pa)
        β = 4.4, # 4.4
        g = [0.0434, 0.0466, 0.2213], # [0.0434, 0.0466, 0.2213]
        ω = 2*π * [10^1, 10^2, 10^3], # 2*π * [10^1, 10^2, 10^3]
        ϵ = 0.9 # 0.9 (AC parameter)
    )
    mesh = Mesh(
        xmin = [-0.3, -0.3, -0.3], # [-0.5, -0.1, -0.1]
        xmax = [ 0.3,  0.3,  0.3], # [ 0.5,  0.1,  0.1]
        tmax = 0.18, # 0.18
        Nx = [100, 100, 100], # [100, 100, 100]
        Nt = 0, # 0
        Co = 0.95, # 0.95 (Courant number)
        flux = "MUSCL", # "LLF", "MUSCL"
        Nf = 6*minimum(length, [mate.g, mate.ω]) + 12 # 6*minimum(length, [mate.g, mate.ω]) + 12
    )
    Sin2(t::Float64, Ω::Float64) = 5e-3 * (sin(Ω*t) - 0.5*sin(2*Ω*t)) .* ((t>=0) - (t>=2*π/Ω))
    src = Source(
        spar = SourcePar(
            type = "Point", # "Plane", "Point"
            dir = 3, # 1
            cmp = 2 # 2
        ),
        Ω = 82.35, # 82.35
        signal = Sin2 # Sin2
    )
    fpara = 9 + src.spar.cmp # field for paraview display
    ana = false # analytical computations
    print("Configuration defined\n")

    # Iterations
    print("Run start (" * Dates.format(now(), "HH:MM") * ")\n")
    @time Run(mate, mesh, src)
    print("Run end (" * Dates.format(now(), "HH:MM") * ")\n")

    # Plot centered slice
    xplot, uplot = Slice(src.spar.dir, mesh)
    hp1 = plot(xplot, uplot[:,10:12], title = string("t = ", mesh.t), xlabel = string("x", src.spar.dir), markershape = :circle, markersize=3)
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
        xth, uth, res = compSol(mate, mesh, src)
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