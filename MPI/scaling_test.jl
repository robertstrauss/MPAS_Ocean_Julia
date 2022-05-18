using MPI
MPI.Init()

comm = MPI.COMM_WORLD
worldrank = MPI.Comm_rank(comm) + 1
root = 1
commsize = MPI.Comm_size(comm)

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
	x = sqrt(n)
	x = Int(round(x))
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

halowidth = 5
ncycles = 5
nsteps = ncycles*halowidth



trials = Int(floor(log2(commsize)))
if worldrank == root
	println("running $trials trials")
end
trialprocs = 2 .^collect(1:trials)
trialmeans = zeros(trials)

for i in 1:trials
	nprocs = trialprocs[i]
	splitcomm = MPI.Comm_split(comm, Int(worldrank>nprocs), worldrank) # only use the necessary ranks for this trial
	if worldrank <= nprocs

		nx, ny = rectangularfactor(nprocs)

		if worldrank == root
			println("dividing $nCellsX x $nCellsX x $(mpasOcean.nVertLevels) ocean among $nprocs processes with $nx x $ny grid")
		end

		cellsInChunk, edgesInChunk, verticesInChunk, cellsFromChunk, cellsToChunk = divide_ocean(fullOcean, halowidth, nx, ny)

		myCells = cellsInChunk[worldrank] # cellsedgesvertices[1]
		myEdges = edgesInChunk[worldrank] # cellsedgesvertices[2]
		myVertices = verticesInChunk[worldrank] # cellsedgesvertices[3]

		mpasOcean = mpas_subset(fullOcean, myCells, myEdges, myVertices)

		kelvinWaveExactNV, kelvinWaveExactSSH, kelvinWaveExactSolution!, boundaryCondition! = kelvinWaveGenerator(mpasOcean, lateralProfilePeriodic)


		# MPI.Barrier(comm)
		# if worldrank == root
		# 	println("kelvin wave initial condition set")
		# end

		if worldrank == root
			println("now simulating for $nsteps steps, communicating every $halowidth steps")
		end

		bench = @benchmark begin
			# simulate
			for f in 1:ncycles
				for h in 1:halowidth
					calculate_normal_velocity_tendency!($mpasOcean)
					update_normal_velocity_by_tendency!($mpasOcean)

					$boundaryCondition!($mpasOcean, exacttime)

					calculate_ssh_tendency!($mpasOcean)
					update_ssh_by_tendency!($mpasOcean)
				end
				exacttime += halowidth*$mpasOcean.dt

				update_halos!($splitcomm, worldrank, $mpasOcean, $cellsFromChunk, $cellsToChunk, $myCells, $myEdges, $myVertices)
			end
		end samples=1 setup=(exacttime=0; $kelvinWaveExactSolution!($mpasOcean, exacttime))

		trialmeans[i] = mean(bench.times)
		if worldrank == root
			println("$(bench.times) ns")
		end

		exactssh = zero(mpasOcean.sshCurrent)
		for iCell in 1:mpasOcean.nCells
			exactssh[iCell,:] = kelvinWaveExactSSH(mpasOcean, iCell, exacttime)
		end

		difference = mpasOcean.sshCurrent .- exactssh
		MaxErrorNorm = norm(difference, Inf)
		L2ErrorNorm = norm(difference/sqrt(float(mpasOcean.nCells)))

		if worldrank == root
			println("error (max, l2): $MaxErrorNorm, $L2ErrorNorm")
		end

	end
	MPI.free(splitcomm)
	MPI.Barrier(comm)
end


if worldrank == root
	pl = Plots.plot(trialprocs, trialmeans, xlabel="number of processors", ylabel="time (ns)",
			title="scaling of distributed simulation of $nCellsX x $nCellsX x $(fullOcean.nVertLevels) ocean kelvin wave for $nsteps steps, MPI every $halowidth",
			titlefontsize=7)
	savedir = "/tmp/rstrauss/mpas-julia"
	mkpath(savedir)
	Plots.savefig(pl, "$savedir/scalingtest_$(Dates.now()).png")
end
