using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
root = 0
commsize = MPI.Comm_size(comm)

include("mode_init/MPAS_Ocean.jl")
include("mode_init/MPAS_OceanHalos.jl")

# if rank == root
println("loading ocean from file...")

fullOcean = MPAS_Ocean("MPAS_O_Shallow_Water/Mesh+Initial_Condition+Registry_Files/Periodic", "base_mesh.nc", "mesh.nc", periodicity="Periodic")

halowidth = 4

# cellsedgesvertices = Array{Int64}(undef, 3)

# recvreq = MPI.Irecv!(cellsedgesvertices, root, rank + commsize, comm)

# if rank == root
cellsInChunk, edgesInChunk, verticesInChunk, cellsFromChunk, cellsToChunk = divide_ocean(fullOcean, halowidth, 2, 2; iChunk = rank+1)



# 	println("scattering ocean pieces to ranks")

# 	sendreqs = []

# 	for destinationrank in 0:commsize-1
# 		iChunk = destinationrank + 1
# 		append!(sendreqs, [MPI.isend([cellsInChunk[iChunk], edgesInChunk[iChunk], verticesInChunk[iChunk]], destinationrank, destinationrank+commsize, comm)])
# 		println("root: sent to rank $destinationrank")
# 	end

# 	MPI.Waitall!([sendreqs...])
# end

# MPI.Wait!(recvreq)

# println("$rank: I recieved!")

myCells = cellsInChunk[rank+1] # cellsedgesvertices[1]
myEdges = edgesInChunk[rank+1] # cellsedgesvertices[2]
myVertices = verticesInChunk[rank+1] # cellsedgesvertices[3]

myOcean = mpas_subset(fullOcean, myCells, myEdges, myVertices)

MPI.Barrier(comm)

println("Ocean distributed between ranks.")

############### Ocean is now distributed between nodes. Let's do some simulation ####################

# set up initial condition

include("mode_init/initial_conditions.jl")

guassianInit!(myOcean)

MPI.Barrier(comm)

println("Initial Condition set.")

# simulate for a while

include("mode_forward/time_steppers.jl")

nFrames = 120

sshOverTime = zeros(mpasOcean.nCells, nFrames)

for f in 1:nFrames
	for h in 1:halowidth
		forward_backward_step!(mpasOcean)
	end
	for (srcchunk, localcells) in cellsFromChunk	
		MPI.Ir
