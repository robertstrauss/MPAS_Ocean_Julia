using MPI
MPI.Init()

comm = MPI.COMM_WORLD
worldrank = MPI.Comm_rank(comm) + 1
root = 1
commsize = MPI.Comm_size(comm)




function rootprint(str)
	if worldrank == root
		println(str)
	end
end


using DelimitedFiles
using LinearAlgebra
using BenchmarkTools
import Plots
import Dates


CODE_ROOT = pwd() * "/"#../"
include(CODE_ROOT * "mode_init/MPAS_Ocean.jl")

include(CODE_ROOT * "mode_init/MPAS_OceanHalos.jl")

include(CODE_ROOT * "mode_init/initial_conditions.jl")
include(CODE_ROOT * "mode_init/exactsolutions.jl")

include(CODE_ROOT * "mode_forward/update_halos.jl")
include(CODE_ROOT * "mode_forward/time_steppers.jl")

# include(CODE_ROOT * "visualization.jl")




function rectangularfactor(n)
	x = Int(round(sqrt(n)))
	if x^2 == n
		return x, x
	end
	while x > 0 && n/x !== round(n/x)
		x -= 1
	end
	return x, Int(round(n/x))
end






nCellsX = 64
fullOcean = MPAS_Ocean(CODE_ROOT * "MPAS_O_Shallow_Water/MPAS_O_Shallow_Water_Mesh_Generation/CoastalKelvinWaveMesh/ConvergenceStudyMeshes",
		       "culled_mesh_$nCellsX.nc", "mesh_$nCellsX.nc", periodicity="NonPeriodic_x", nvlevels=100)

lYedge = maximum(fullOcean.yEdge) - minimum(fullOcean.yEdge)
lateralProfilePeriodic(y) = 1e-3*cos(y/lYedge * 4 * pi)


rootprint("preparing kelvin wave")
kelvinWaveExactNV, kelvinWaveExactSSH, kelvinWaveExactSolution!, boundaryCondition! = kelvinWaveGenerator(fullOcean, lateralProfilePeriodic)


halowidth = 5
ncycles = 2
nsteps = ncycles*halowidth


ntests = Int(floor(log2(commsize))) # number of levels (counts of processors) to run for
nsamples = 3 # times to rerun sim for each level
proccounts = 2 .^collect(1:ntests)
rootprint("running tests for $prouccounts counts of processors")
wctime = zeros((ntests,nsamples))


partitiondir = CODE_ROOT * "scaling_test/graphparts/"

function runtrials()
	for itest in 1:ntests
		nprocs = proccounts[itest]
		rootprint("splitting com")
		splitcomm = MPI.Comm_split(comm, Int(worldrank>nprocs), worldrank) # only use the necessary ranks for this trial
		if worldrank <= nprocs

			rootprint("doing rect facto")
			nx, ny = rectangularfactor(nprocs)

			# partitionfile = partitiondir * "graph.info.part.$nprocs"
			# rootprint("partitioning $nCellsX x $nCellsX x $(fullOcean.nVertLevels) ocean among $nprocs processes with $(partitionfile)")
			# partitions = dropdims(readdlm(partitionfile), dims=2)
			# mycells = findall(proc -> proc == worldrank-1, partitions)

			cellsInChunk, edgesInChunk, verticesInChunk, cellsFromChunk, cellsToChunk = divide_ocean(fullOcean, halowidth, nx, ny)


			myCells = cellsInChunk[worldrank] # cellsedgesvertices[1]
			myEdges = edgesInChunk[worldrank] # cellsedgesvertices[2]
			myVertices = verticesInChunk[worldrank] # cellsedgesvertices[3]

			mpasOcean = mpas_subset(fullOcean, myCells, myEdges, myVertices)



			# MPI.Barrier(comm)
			# if worldrank == root
			# 	println("kelvin wave initial condition set")
			# end

			rootprint("now simulating for $nsteps steps, communicating every $halowidth steps")

			exacttime=0; kelvinWaveExactSolution!(mpasOcean, exacttime)
			for jsample in 1:nsamples
				sampletime = @elapsed begin
					# simulate
					for f in 1:ncycles
						for h in 1:halowidth
							calculate_normal_velocity_tendency!(mpasOcean)
							update_normal_velocity_by_tendency!(mpasOcean)

							boundaryCondition!(mpasOcean, exacttime)

							calculate_ssh_tendency!(mpasOcean)
							update_ssh_by_tendency!(mpasOcean)
						end
						exacttime += halowidth*mpasOcean.dt

						update_halos!(splitcomm, worldrank, mpasOcean, cellsFromChunk, cellsToChunk, myCells, myEdges, myVertices)
					end
				end
				wctime[itest, jsample] = sampletime
			end # samples=1 setup=(exacttime=0; $kelvinWaveExactSolution!($mpasOcean, exacttime))

			rootprint("$(wctime[itest,:]) ns")

			rootprint("setting up exact sol")
			exactssh = zero(mpasOcean.sshCurrent)
			for iCell in 1:mpasOcean.nCells
				exactssh[iCell,:] .= kelvinWaveExactSSH(mpasOcean, iCell, nsteps*mpasOcean.dt)
			end

			rootprint("calculating error")
			difference = mpasOcean.sshCurrent .- exactssh
			MaxErrorNorm = norm(difference, Inf)
			L2ErrorNorm = norm(difference/sqrt(float(mpasOcean.nCells)))

			rootprint("error (max, l2): $MaxErrorNorm, $L2ErrorNorm")

		end
		rootprint("freeing splitcomm")
		MPI.free(splitcomm)
		MPI.Barrier(comm)
	end


	if worldrank == root

		fpath = "/global/cscratch1/sd/rstrauss/scaling_test_1/kelvinwave/resolution$(nCellsX)x$(nCellsX)/steps$(nsteps)"
		fname = fpath * "/$(Dates.now()).txt"
		mkpath(fpath)
		open(fname, "w") do io
			writedlm(io, [trialprocs, trialmeans], ',')
		end
	#=	pl = Plots.plot(trialprocs, trialmeans, xlabel="number of processors", ylabel="time (ns)",
				title="scaling of distributed simulation of $nCellsX x $nCellsX x $(fullOcean.nVertLevels) ocean kelvin wave for $nsteps steps, MPI every $halowidth",
				titlefontsize=7)
		savedir = "/tmp/rstrauss/mpas-julia"
		mkpath(savedir)
		Plots.savefig(pl, "$savedir/scalingtest_$(Dates.now()).png")
		=#
	end
end
