using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm) + 1
root = 1
commsize = MPI.Comm_size(comm)

using DelimitedFiles
using LinearAlgebra




CODE_ROOT = pwd() * "/"#../"
include(CODE_ROOT * "mode_init/MPAS_Ocean.jl")


println("loading ocean from file...")

cd(CODE_ROOT)
# fullOcean = MPAS_Ocean("MPAS_O_Shallow_Water/Mesh+Initial_Condition+Registry_Files/Periodic", "base_mesh.nc", "mesh.nc", periodicity="Periodic")
nCellsX = 64
fullOcean = MPAS_Ocean(CODE_ROOT * "MPAS_O_Shallow_Water/MPAS_O_Shallow_Water_Mesh_Generation/CoastalKelvinWaveMesh/ConvergenceStudyMeshes",
		       "culled_mesh_$nCellsX.nc", "mesh_$nCellsX.nc", periodicity="NonPeriodic_x")
halowidth = 5



include(CODE_ROOT * "mode_init/MPAS_OceanHalos.jl")

cellsInChunk, edgesInChunk, verticesInChunk, cellsFromChunk, cellsToChunk = divide_ocean(fullOcean, halowidth, 2, 2)#; iChunk = rank)

myCells = cellsInChunk[rank] # cellsedgesvertices[1]
myEdges = edgesInChunk[rank] # cellsedgesvertices[2]
myVertices = verticesInChunk[rank] # cellsedgesvertices[3]

mpasOcean = mpas_subset(fullOcean, myCells, myEdges, myVertices)

MPI.Barrier(comm)

if rank == root
	println("ocean distributed between ranks.")
end



include(CODE_ROOT * "mode_init/initial_conditions.jl")
include(CODE_ROOT * "mode_init/exactsolutions.jl")


lYedge = maximum(fullOcean.yEdge) - minimum(fullOcean.yEdge)
lateralProfilePeriodic(y) = 1e-3*cos(y/lYedge * 4 * pi)
kelvinWaveExactNV, kelvinWaveExactSSH, kelvinWaveExactSolution!, boundaryCondition! = kelvinWaveGenerator(mpasOcean, lateralProfilePeriodic)


include(CODE_ROOT * "mode_forward/update_halos.jl")
include(CODE_ROOT * "mode_forward/time_steppers.jl")

include(CODE_ROOT * "visualization.jl")

function kelvinwave(n, writeoutput=false, displayplots=false)
	# sshOverTime = zeros(n+1, mpasOcean.nCells)
	# sshOverTime[1,:] = mpasOcean.sshCurrent
	exacttime = 0
	kelvinWaveExactSolution!(mpasOcean, exacttime)

	MPI.Barrier(comm)
	if rank == root
		println("kelvin wave initial condition set")
	end

	# simulate for a while

	for f in 1:n
		for h in 1:halowidth
			calculate_normal_velocity_tendency!(mpasOcean)
			update_normal_velocity_by_tendency!(mpasOcean)

			boundaryCondition!(mpasOcean, exacttime)

			calculate_ssh_tendency!(mpasOcean)
			update_ssh_by_tendency!(mpasOcean)
		end
		exacttime += halowidth*mpasOcean.dt

		update_halos!(comm, rank, mpasOcean, cellsFromChunk, cellsToChunk, myCells, myEdges, myVertices)

		# sshOverTime[f+1,:] = mpasOcean.sshCurrent
		#
		# if rank == root
		# 	println("$rank: iteration $f of $n complete")
		# end
	end

	exactssh = zero(mpasOcean.sshCurrent)
	for iCell in 1:mpasOcean.nCells
		exactssh[iCell] = kelvinWaveExactSSH(mpasOcean, iCell, exacttime)
	end

	difference = mpasOcean.sshCurrent .- exactssh
	MaxErrorNorm = norm(difference, Inf)
	L2ErrorNorm = norm(difference/sqrt(float(mpasOcean.nCells)))

	if displayplots# && rank == root
		fig, axes = subplots(1,3, figsize=(9,3))
		_, ax, _, _ = heatMapMesh(mpasOcean, mpasOcean.sshCurrent, fig=fig, ax=axes[1])
        axes[1].title.set_text("numerical solution")

		_, ax, _, _ = heatMapMesh(mpasOcean, exactssh, fig=fig, ax=axes[2])
        axes[2].title.set_text("exact solution")

        _, ax, _, _ = heatMapMesh(mpasOcean, difference, fig=fig, ax=axes[3])
        axes[3].title.set_text("difference")

		fig.suptitle("rank $rank")
        display(fig)
		fig.savefig("/tmp/kelvinwave/ssh$rank.png")
	end


	return mpasOcean.gridSpacingMagnitude, MaxErrorNorm, L2ErrorNorm
end



kelvinwave(10, false, true)
