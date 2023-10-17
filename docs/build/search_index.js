var documenterSearchIndex = {"docs":
[{"location":"man.html","page":"Manual","title":"Manual","text":"The reference documentation.","category":"page"},{"location":"man.html#Getting-started","page":"Manual","title":"Getting started","text":"","category":"section"},{"location":"man.html#Running-a-simulation","page":"Manual","title":"Running a simulation","text":"","category":"section"},{"location":"man.html#Running-examples","page":"Manual","title":"Running examples","text":"","category":"section"},{"location":"api.html","page":"API","title":"API","text":"The application programming interface.","category":"page"},{"location":"api.html#Functions","page":"API","title":"Functions","text":"","category":"section"},{"location":"api.html","page":"API","title":"API","text":"CurrentModule = SoftSol3D\nDocTestSetup = quote\n    using SoftSol3D\nend","category":"page"},{"location":"api.html","page":"API","title":"API","text":"Modules = [SoftSol3D, SoftSol3D.ConfigPar]\nPrivate = false\nOrder   = [:module, :function, :type]","category":"page"},{"location":"api.html#SoftSol3D.LocalMaxSpeeds-Tuple{Vector{Float64}, Material, Vector{Float64}}","page":"API","title":"SoftSol3D.LocalMaxSpeeds","text":"LocalMaxSpeeds(Q::Array{Float64,1}, mate::Material, cmaxloc::Array{Float64,1})\n\nReturns the local maximum characteristic speeds along x, y, z at the state vector Q with the parameters of the argument.\n\n\n\n\n\n","category":"method"},{"location":"api.html#SoftSol3D.Run-Tuple{Material, Mesh, Source}","page":"API","title":"SoftSol3D.Run","text":"Run(mate::Material, mesh::Mesh, src::Source)\n\nExecutes the main loop for the parameters specified in the argument. The numerical solution  at the final time data.t is stored in data.Qp. This data vector of size 6N+13 has the following structure:\n\nbf Q = left(F_11 F_12 F_13 dots  F_33 v_1 v_2 v_3 S_1^textv_11 dots  S_N^textv_33right)^top\n\nwhere F_ij are the deformation gradient components, v_i are the velocity components, and S_ell^textv_ij are the memory variables.\n\n\n\n\n\n","category":"method"},{"location":"api.html#SoftSol3D.Slice-Tuple{Int64, Mesh}","page":"API","title":"SoftSol3D.Slice","text":"x, u = Slice(dir::Int64, mesh::Mesh, data::Data)\n\nReturns data data.Qp sliced along direction dir in output u. Spatial coordinates mesh.X in output x undergo the same operation.\n\n\n\n\n\n","category":"method"},{"location":"api.html#SoftSol3D.compSol-Tuple{Material, Mesh, Source}","page":"API","title":"SoftSol3D.compSol","text":"xth, uth, res = compSol(mate::Material, mesh::Mesh, src::Source)\n\nComputes analytical solutions via FFT.\n\n\n\n\n\n","category":"method"},{"location":"api.html#SoftSol3D.ConfigPar","page":"API","title":"SoftSol3D.ConfigPar","text":"Defines the data types using Parameters.jl.\n\n\n\n\n\n","category":"module"},{"location":"api.html#SoftSol3D.ConfigPar.initSource-Tuple{Material, Mesh, Source}","page":"API","title":"SoftSol3D.ConfigPar.initSource","text":"pts, coeffs = initSource(mate::Material, mesh::Mesh, src::Source)\n\nSource configuration. Returns source points pts and corresponding body force coefficients coeffs.\n\n\n\n\n\n","category":"method"},{"location":"api.html#SoftSol3D.ConfigPar.Material","page":"API","title":"SoftSol3D.ConfigPar.Material","text":"Material parameters.\n\n\n\n\n\n","category":"type"},{"location":"api.html#SoftSol3D.ConfigPar.Mesh","page":"API","title":"SoftSol3D.ConfigPar.Mesh","text":"Computational parameters.\n\n\n\n\n\n","category":"type"},{"location":"api.html#SoftSol3D.ConfigPar.Source","page":"API","title":"SoftSol3D.ConfigPar.Source","text":"Source parameters (loading).\n\n\n\n\n\n","category":"type"},{"location":"api.html#SoftSol3D.ConfigPar.SourcePar","page":"API","title":"SoftSol3D.ConfigPar.SourcePar","text":"Source sub-parameters.\n\n\n\n\n\n","category":"type"},{"location":"api.html#Index","page":"API","title":"Index","text":"","category":"section"},{"location":"api.html","page":"API","title":"API","text":"","category":"page"},{"location":"index.html#SoftSol3D.jl","page":"Home","title":"SoftSol3D.jl","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"A Julia package for the 3D incompressible quasi-linear viscoelasticity equations.","category":"page"},{"location":"index.html#Version","page":"Home","title":"Version","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"using Dates # hide\nprint(\"Built UTC \", now()) # hide","category":"page"},{"location":"index.html#Author","page":"Home","title":"Author","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"Harold Berjamin, University of Galway","category":"page"},{"location":"index.html#License","page":"Home","title":"License","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"SoftSol3D is freely distributed under the MIT license.","category":"page"},{"location":"index.html#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"SoftSol3D.jl is an unregistered package, and is simply installed by downloading and extracting the files. Then, open Julia's REPL, type ] and run","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"pkg> add /home/JuliaUser/Downloads/SoftSol3D","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"where /home/JuliaUser/Downloads/SoftSol3D is the path to the SoftSol3D folder.","category":"page"},{"location":"index.html#Citation","page":"Home","title":"Citation","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"H. Berjamin (2023). \"Computation of viscoelastic shear shock waves using finite volume schemes with artificial compressibility\", Arxiv preprint 2310.04355. doi:10.48550/arXiv.2310.04355","category":"page"},{"location":"index.html#Acknowledgments","page":"Home","title":"Acknowledgments","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"This project began during a Postdoc project at the University of Galway. It has received funding from the European Union's Horizon 2020 research and innovation programme under Grant No. TBI-WAVES—H2020-MSCA-IF-2020 Project No. 101023950.","category":"page"}]
}
