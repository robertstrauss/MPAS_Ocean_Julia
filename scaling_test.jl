using MPI
MPI.Init()

comm = MPI.COMM_WORLD
worldrank = MPI.Comm_rank(comm) + 1
root = 1
commsize = MPI.Comm_size(comm)


using DelimitedFiles
using DataFrames
using CSV
using LinearAlgebra
using BenchmarkTools
# import Plots
import Dates

using TimerOutputs


function rootprint(str)
    if worldrank == root
        println("$(Dates.Time(Dates.now())): ", str)
    end
end



CODE_ROOT = pwd() * "/"
include(CODE_ROOT * "mode_init/MPAS_Ocean.jl")
include(CODE_ROOT * "mode_init/MPAS_OceanHalos.jl")
include(CODE_ROOT * "mode_init/exactsolutions.jl")

include(CODE_ROOT * "mode_forward/update_halos.jl")
include(CODE_ROOT * "mode_forward/time_steppers.jl")




function globalToLocal(globalcells, mycells)
    localcells = findall(iCell -> iCell in globalcells, mycells)
    return localcells
end


function runtests(proccounts, fname, partitiondir; nsamples=6, nCellsX=64, halowidth=5, ncycles=2, nvlevels=100)

    meshpath = CODE_ROOT * "MPAS_O_Shallow_Water/ConvergenceStudyMeshes/CoastalKelvinWave" #MPAS_O_Shallow_Water_Mesh_Generation/CoastalKelvinWaveMesh/ConvergenceStudyMeshes"

    # kelvin wave initial condition and exact solution functions
    if worldrank == root
        mesh_file = NCDatasets.Dataset("$meshpath/mesh_$(nCellsX)x$(nCellsX).nc", "r", format=:netcdf4)
        # fullOcean = MPAS_Ocean(meshpath, "culled_mesh_$nCellsX.nc", "mesh_$nCellsX.nc", periodicity="NonPeriodic_x", nvlevels=nvlevels)
        gravity = 9.8 # fullOcean.gravity
        # dt = fullOcean.dt
        dcEdge = mesh_file["dcEdge"]
        yEdge = mesh_file["yEdge"]
        lYedge = maximum(yEdge) - minimum(yEdge)
        # lYedge = maximum(fullOcean.yEdge) - minimum(fullOcean.yEdge)
        # meanCoriolisParameterf = sum(fullOcean.fEdge) / length(fullOcean.fEdge)
        # meanFluidThicknessH = sum(fullOcean.bottomDepth)/length(fullOcean.bottomDepth)
        meanCoriolisParameterf = 0.0001 # sum(fEdge) / length(fEdge)
        meanFluidThicknessH = 1000.0 # sum(bottomDepth) / length(bottomDepth)
        c = sqrt(gravity*meanFluidThicknessH)
        rossbyRadiusR = c/meanCoriolisParameterf

        courantNumber = 0.2
        dt = courantNumber * minimum(dcEdge) / c
    else
        lYedge = 0.0
        meanCoriolisParameterf = 0.0
        meanFluidThicknessH = 0.0
        c = 0.0
        rossbyRadiusR = 0.0
        gravity = 0.0
        dt = 0.0
    end
    MPI.Barrier(comm)
    lYedge, meanCoriolisParameterf, meanFluidThicknessH, c, rossbyRadiusR, gravity, dt = MPI.bcast([lYedge, meanCoriolisParameterf, meanFluidThicknessH, c, rossbyRadiusR, gravity, dt], root-1, comm)

    lateralProfilePeriodic(y) = 1e-6*cos(y/lYedge * 4 * pi)
    lateralProfile = lateralProfilePeriodic

    # kelvinWaveExactNormalVelocity, kelvinWaveExactSSH, kelvinWaveExactSolution!, boundaryCondition! = kelvinWaveGenerator(mpasOcean, lateralProfile)
    # """

    function kelvinWaveExactNormalVelocity(mpasOcean, iEdge, t=0)
        v = c * lateralProfile.(mpasOcean.yEdge[iEdge] .+ c*t) .* exp.(-mpasOcean.xEdge[iEdge]/rossbyRadiusR)
        return v .* sin.(mpasOcean.angleEdge[iEdge])
        end

    function kelvinWaveExactSSH(mpasOcean, iCell, t=0)
        return - meanFluidThicknessH * lateralProfile.(mpasOcean.yCell[iCell] .+ c*t) .* exp.(-mpasOcean.xCell[iCell]/rossbyRadiusR)
    end

    function kelvinWaveExactSolution!(mpasOcean, t=0)
        # just calculate exact solution once then copy it to lower layers
        sshperlayer = kelvinWaveExactSSH(mpasOcean, collect(1:mpasOcean.nCells), t) / mpasOcean.nVertLevels
        for k in 1:mpasOcean.nVertLevels
            mpasOcean.layerThickness[k,:] .+= sshperlayer # add, don't replace, thickness contributes to level depth
        end

        nvperlayer = kelvinWaveExactNormalVelocity(mpasOcean, collect(1:mpasOcean.nEdges), t) / mpasOcean.nVertLevels
        for k in 1:mpasOcean.nVertLevels
            mpasOcean.normalVelocityCurrent[k,:] .= nvperlayer
        end
    end

    function boundaryCondition!(mpasOcean, t)
        for iEdge in 1:mpasOcean.nEdges
            if mpasOcean.boundaryEdge[iEdge] == 1.0
                mpasOcean.normalVelocityCurrent[:,iEdge] .= kelvinWaveExactNormalVelocity(mpasOcean, iEdge, t)/mpasOcean.nVertLevels
            end
        end
    end
    # """


    rootprint("running tests for $proccounts counts of processors with $nsamples samples per test of $(halowidth*ncycles) steps each")
    rootprint("halos: $halowidth wide")


    df = DataFrame([Int64[], collect([Float64[] for i in 1:nsamples])...,     collect([Float64[] for i in 1:nsamples])...,     Float64[], Float64[]],
               ["procs", collect(["sim_time$i" for i in 1:nsamples])...,  collect(["mpi_time$i" for i in 1:nsamples])...,  "max_error", "l2_error"])


    for nprocs in proccounts
        splitcomm = MPI.Comm_split(comm, Int(worldrank>nprocs), worldrank) # only use the necessary ranks for this trial
        if worldrank <= nprocs

	    myoceancells = "All"
            if nprocs > 1
		partitionfile = partitiondir * "graph.info.part.$nprocs"
		rootprint("partitioning $nCellsX x $nCellsX x $(nvlevels) ocean among $nprocs processes with $(partitionfile)")
		partitions = dropdims(readdlm(partitionfile), dims=2)

		cellsOnCell = Array{Int64,2}(replace(readdlm(partitiondir * "graph.info")[2:end,:], ""=>0))

		haloMesh = MPAS_OceanHalos(splitcomm, cellsOnCell, partitions)

                rootprint("cell count $(length(haloMesh.myCells))")

		myoceancells = haloMesh.myCells
            end

	    mpasOcean = MPAS_Ocean(meshpath, "culled_mesh_$(nCellsX)x$(nCellsX).nc", "mesh_$(nCellsX)x$(nCellsX).nc", periodicity="NonPeriodic_x", nvlevels=nvlevels, cells=myoceancells)
            mpasOcean.dt = dt #1e-2 * minimum(mpasOcean.dcEdge) / c
            rootprint("dt: $(mpasOcean.dt)")


            if nprocs > 1
		fillEdgeArraysFromCells!(haloMesh, mpasOcean.edgesOnCell)
		generateBufferArrays!(haloMesh, mpasOcean)

		rootprint("rank $(worldrank): cellsToMyRankFrom:  $(collect(["$k -> $(length(lc))" for (k, lc) in haloMesh.cellsToMyRankFrom]))")
		# println("rank $(worldrank): haloBuffersCells: $(collect(["$k -> $(length(lc))" for (k, lc) in haloMesh.haloBuffersCells]))")
	    end

            rootprint("now simulating for $(halowidth*ncycles*nsamples) steps, communicating every $halowidth steps, timing $(nsamples) samples of $(halowidth*ncycles) steps")

            sampletimessim = zeros(nsamples)
            sampletimesmpi = zeros(nsamples)
            exacttime=0; kelvinWaveExactSolution!(mpasOcean, exacttime)

            for jsample in 1:nsamples
                for f in 1:ncycles
                    # simulate
                    sampletimessim[jsample] += @elapsed begin
                        for h in 1:halowidth
                            calculate_normal_velocity_tendency!(mpasOcean)
                            update_normal_velocity_by_tendency!(mpasOcean)

                            exacttime += mpasOcean.dt
                            boundaryCondition!(mpasOcean, exacttime)

                            calculate_thickness_tendency!(mpasOcean)
                            update_thickness_by_tendency!(mpasOcean)
                        end
                    end
                    if nprocs > 1
                        sampletimesmpi[jsample] += @elapsed begin
			    update_halos!(splitcomm, mpasOcean, haloMesh) 
                        end
                    end
                end
		Base.GC.gc()
            end

            rootprint("generating exact sol")
            exactssh = zero(mpasOcean.sshCurrent)
            for iCell in 1:mpasOcean.nCells
                exactssh[iCell] = kelvinWaveExactSSH(mpasOcean, iCell, halowidth*ncycles*mpasOcean.dt)
            end

            rootprint("calculating error")
            mpasOcean.sshCurrent = dropdims(sum(mpasOcean.layerThickness, dims=1), dims=1) - mpasOcean.bottomDepth
            difference = mpasOcean.sshCurrent .- exactssh
            MaxErrorNorm = norm(difference, Inf)
            L2ErrorNorm = norm(difference/sqrt(float(mpasOcean.nCells)))

            push!(df, [nprocs, sampletimessim..., sampletimesmpi..., MaxErrorNorm, L2ErrorNorm])

            rootprint(df[df.procs .<= nprocs, :])
            if worldrank == root
                CSV.write(fname, df)
            end
            rootprint("saved at $fname")

        end
        MPI.free(splitcomm)
        MPI.Barrier(comm)

    end



    return
end

halowidth = 1
ncycles = 10

if length(ARGS) < 2
    println("Invalid call signature! Format: ./scaling_test.jl [nCellsX] [nSamples] [list of n procs (optional)]")
end
nCellsX = parse(Int64, ARGS[1])
nsamples = parse(Int64, ARGS[2])
if length(ARGS) > 2
    proccounts = parse.(Int64, split(ARGS[3], ','))
else
    proccounts = Int.(round.( 2 .^collect(0:log2(commsize)) ))
end
nvlevels = 100

partitiondir = CODE_ROOT * "scaling_test/graphparts/$(nCellsX)x$(nCellsX)/"
fpath = CODE_ROOT * "output/kelvinwave/resolution$(nCellsX)x$(nCellsX)/procs$(commsize)/steps$(halowidth*ncycles)/nvlevels$(nvlevels)"
mkpath(fpath)
fname = fpath * "/$(Dates.now()).txt"
rootprint("results to: $fname")

runtests(proccounts, fname, partitiondir; nCellsX=nCellsX, nsamples=nsamples, halowidth=halowidth, ncycles=ncycles, nvlevels=nvlevels)
