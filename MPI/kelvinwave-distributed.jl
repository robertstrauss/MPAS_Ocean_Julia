using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm) + 1
root = 1
commsize = MPI.Comm_size(comm)


using DelimitedFiles

CODE_ROOT = pwd() * "/"#../"

include(CODE_ROOT * "mode_init/MPAS_Ocean.jl")
include(CODE_ROOT * "mode_init/MPAS_OceanHalos.jl")


include(CODE_ROOT * "mode_init/initial_conditions.jl")

include(CODE_ROOT * "mode_forward/time_steppers.jl")

include(CODE_ROOT * "mode_init/exactsolutions.jl")


println("loading ocean from file...")

cd(CODE_ROOT)


F = Float64
I = Int64


# fullOcean = MPAS_Ocean("MPAS_O_Shallow_Water/Mesh+Initial_Condition+Registry_Files/Periodic", "base_mesh.nc", "mesh.nc", periodicity="Periodic")
nCellsX = 64::I
fullOcean = MPAS_Ocean(CODE_ROOT * "MPAS_O_Shallow_Water/MPAS_O_Shallow_Water_Mesh_Generation/CoastalKelvinWaveMesh/ConvergenceStudyMeshes",
		       "culled_mesh_$nCellsX.nc", "mesh_$nCellsX.nc", periodicity="NonPeriodic_x")
halowidth = 5::I

cellsInChunk, edgesInChunk, verticesInChunk, cellsFromChunk, cellsToChunk = divide_ocean(fullOcean, halowidth, 2, 2)#; iChunk = rank)

myCells = cellsInChunk[rank] # cellsedgesvertices[1]
myEdges = edgesInChunk[rank] # cellsedgesvertices[2]
myVertices = verticesInChunk[rank] # cellsedgesvertices[3]

mpasOcean = mpas_subset(fullOcean, myCells, myEdges, myVertices)

MPI.Barrier(comm)

if rank == root
	println("ocean distributed between ranks.")
end




lYedge = maximum(fullOcean.yEdge) - minimum(fullOcean.yEdge)
lateralProfilePeriodic(y) = 1e-3*cos(y/lYedge * 4 * pi)
kelvinWaveExactNV, kelvinWaveExactSSH, kelvinWaveExactSolution!, boundaryCondition! = kelvinWaveGenerator(mpasOcean, lateralProfilePeriodic)


exacttime = 0.0::F


function stepAndSync()
	global exacttime
	# simulate until the halo areas are all invalid and need to be updated
	for h in 1:halowidth
		# forward_backward_step!(mpasOcean)
		calculate_normal_velocity_tendency!(mpasOcean)
		update_normal_velocity_by_tendency!(mpasOcean)

		boundaryCondition!(mpasOcean, exacttime)

		calculate_ssh_tendency!(mpasOcean)
		update_ssh_by_tendency!(mpasOcean)

		exacttime += mpasOcean.dt
	end

	### request cells in my halo from chunks with those cells
	halobufferssh = [] # temporarily stores new halo ssh
	halobuffernv = [] # temporarily stores new halo normal velocity
	recreqs = []
	for (srcchunk, localcells) in cellsFromChunk[rank]
		newhalossh = Array{eltype(mpasOcean.sshCurrent)}(undef, length(localcells))
		append!(halobufferssh, [newhalossh])
		reqssh = MPI.Irecv!(newhalossh, srcchunk-1, 0, comm) # tag 0 for ssh
		append!(recreqs, [reqssh])

		localedges = collect(Set(mpasOcean.edgesOnCell[:,localcells]))
		newhalonv = Array{eltype(mpasOcean.normalVelocityCurrent)}(undef, length(localedges))
		append!(halobuffernv, [newhalonv])
		reqnv = MPI.Irecv!(newhalonv, srcchunk-1, 1, comm) # tag 1 for norm vel
		append!(recreqs, [reqnv])
	end

	MPI.Barrier(comm)
	### send cells in main non-halo area to chunks that need them for their halo
	sendreqs = []
	for (dstchunk, localcells) in cellsToChunk[rank]
		reqssh = MPI.Isend(mpasOcean.sshCurrent[localcells], dstchunk-1, 0, comm)
		append!(sendreqs, [reqssh])

		localedges = collect(Set(mpasOcean.edgesOnCell[:,localcells])) # Set to remove duplicates
		order = sortperm(myEdges[localedges])
		reqnv = MPI.Isend(mpasOcean.normalVelocityCurrent[localedges[order]], dstchunk-1, 1, comm)
		append!(sendreqs, [reqnv])
	end

	### copy the recieved data into the ocean's halo
	if rank == root
	# 	println("halo buffer before: ", halobuffernv[1][1:10])
	end
	MPI.Barrier(comm)
	MPI.Waitall!([recreqs..., sendreqs...])
	if rank == root
	# 	println("halo buffer after: ", halobuffernv[1][1:10])
	end
	MPI.Barrier(comm)
	for (i, (_, localcells)) in enumerate(cellsFromChunk[rank])
		mpasOcean.sshCurrent[localcells] = halobufferssh[i]
		localedges = collect(Set(mpasOcean.edgesOnCell[:,localcells]))
		order = sortperm(myEdges[localedges])
		mpasOcean.normalVelocityCurrent[localedges[order]] = halobuffernv[i]
	end
end


function kelvinwave(n, writeoutput=false)
	sshOverTime = zeros(n+1, mpasOcean.nCells)
	sshOverTime[1,:] = mpasOcean.sshCurrent


# 	exacttime = 0
	kelvinWaveExactSolution!(mpasOcean, exacttime)

	MPI.Barrier(comm)

# 	if rank == root
		println("$rank: initial condition set.")
	# end

	# simulate for a while

	for f in 1:n
		stepAndSync()

		# exacttime += halowidth*dt

		sshOverTime[f+1,:] = mpasOcean.sshCurrent

		if rank == root
			println("$rank: iteration $f of $n complete")
		end
	end

	print("$rank: total error (ssh): $(sum([abs(mpasOcean.sshCurrent[iCell] - kelvinWaveExactSSH(mpasOcean, iCell, exacttime)) for iCell in 1:mpasOcean.nCells]))")

	if writeoutput
		if rank == root
			println("simulation complete, writing output arrays")
		end


		open("out/dist_test_1_ssh_rank_$rank.txt", "w") do io
			writedlm(io, sshOverTime)
		end

		println("$rank: file written")
	end
end

# # force compilation
stepAndSync()
kelvinwave(1)

using Profile

@profile stepAndSync()

if rank == 1
	print("============= step and sync profile =================")
	Profile.print()
end
