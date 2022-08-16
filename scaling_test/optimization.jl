
CODE_ROOT = pwd() * "/"#../"

include(CODE_ROOT * "mode_init/MPAS_Ocean.jl")
include(CODE_ROOT * "mode_init/MPAS_OceanHalos.jl")

include(CODE_ROOT * "mode_init/initial_conditions.jl")
include(CODE_ROOT * "mode_forward/time_steppers.jl")
include(CODE_ROOT * "mode_init/exactsolutions.jl")

using BenchmarkTools

using DelimitedFiles

rank = 1
worldrank = 1
root = 1

function rootprint(text)
	if rank == root
		println(text)
	end
end


nprocs = 4

# fullOcean = MPAS_Ocean("MPAS_O_Shallow_Water/Mesh+Initial_Condition+Registry_Files/Periodic", "base_mesh.nc", "mesh.nc", periodicity="Periodic")
nCellsX = 64
println("loading ocean from file")
fullOcean = MPAS_Ocean(
				CODE_ROOT * "MPAS_O_Shallow_Water/MPAS_O_Shallow_Water_Mesh_Generation/CoastalKelvinWaveMesh/ConvergenceStudyMeshes",
		       "culled_mesh_$nCellsX.nc", "mesh_$nCellsX.nc", nvlevels=1, periodicity="NonPeriodic_x")
halowidth = 5


partitiondir = CODE_ROOT * "scaling_test/graphparts/$(nCellsX)x$(nCellsX)/"
# fpath = CODE_ROOT * "output/kelvinwave/resolution$(nCellsX)x$(nCellsX)/steps$(halowidth*ncycles)"


println("dividing ocean")
# cellsInChunk, edgesInChunk, verticesInChunk, cellsFromChunk, cellsToChunk = divide_ocean(fullOcean, halowidth, 2, 2)#; iChunk = rank)
function partitionorig(partitionfile::String = partitiondir * "graph.info.part.$nprocs")
	# rootprint("partitioning $nCellsX x $nCellsX x $(fullOcean.nVertLevels) ocean among $nprocs processes with $(partitionfile)")
	partitions = dropdims(readdlm(partitionfile), dims=2)
	mycells = findall(proc -> proc == worldrank-1, partitions)

	# rootprint("halos: $halowidth wide")
	haloCells = grow_halo(fullOcean, mycells, halowidth)
	# rootprint("\t halo $haloCells")

	cellsFromChunk = Dict{Int64, Array}()
	for iCell in haloCells
		chunk = partitions[iCell] + 1
		if chunk in keys(cellsFromChunk)
			push!(cellsFromChunk[chunk], iCell)
		else
			cellsFromChunk[chunk] = [iCell]
		end
	end

	# tell neighbors which cells they should send me
	cellsToChunk = Dict{Int64, Array}()
	tag = 2
	# sendreqs = []; recvreqs = []
	for (chunk, cells) in cellsFromChunk
		# push!(sendreqs, MPI.Isend(cells, chunk-1, tag, splitcomm))
		#
		# numcells = MPI.Get_count(MPI.Probe(chunk-1, tag, splitcomm), Int64)
		# cellsToChunk[chunk] = Array{Int64,1}(undef, numcells)
		# push!(recvreqs, MPI.Irecv!(cellsToChunk[chunk], chunk-1, tag, splitcomm))
	end
	# MPI.Waitall!([sendreqs..., recvreqs...])
end


function partition(partitionfile::String = partitiondir * "graph.info.part.$nprocs")
	# rootprint("partitioning $nCellsX x $nCellsX x $(fullOcean.nVertLevels) ocean among $nprocs processes with $(partitionfile)")
	partitions::Array{Int64,1} = dropdims((readdlm(partitionfile)), dims=2)
	# partitions = Array{Int64,1}()
	# println("arr size ", sizeof(partitions))
	mycells = findall(proc -> proc == worldrank-1, partitions)

	# rootprint("halos: $halowidth wide")
	haloCells = grow_halo(fullOcean, mycells, halowidth)
	# haloCells = Array{Int64, 1}()
	# rootprint("\t halo $haloCells")

	cellsFromChunk = Dict{Int64, Array{Int64,1}}()
	for iCell in haloCells
		chunk = partitions[iCell] + 1
		if chunk in keys(cellsFromChunk)
			push!(cellsFromChunk[chunk], iCell)
		else
			cellsFromChunk[chunk] = [iCell]
		end
	end

	# tell neighbors which cells they should send me
	cellsToChunk = Dict{Int64, Array{Int64,1}}()
	tag = 2
	# sendreqs = []; recvreqs = []
	for (chunk::Int64, cells::Array{Int64}) in cellsFromChunk
		# push!(sendreqs, MPI.Isend(cells, chunk-1, tag, splitcomm))

		# numcells = MPI.Get_count(MPI.Probe(chunk-1, tag, splitcomm), Int64)
		# cellsToChunk[chunk] = Array{Int64,1}(undef, numcells)
		# push!(recvreqs, MPI.Irecv!(cellsToChunk[chunk], chunk-1, tag, splitcomm))
	end
	# MPI.Waitall!([sendreqs..., recvreqs...])
end
# partition(partitiondir * "graph.info.part.$nprocs")

@code_warntype partition(partitiondir * "graph.info.part.$nprocs")

@btime partition(partitiondir * "graph.info.part.$nprocs")

@btime partitionorig()


include(CODE_ROOT * "mode_init/MPAS_OceanHalos.jl")

@code_warntype grow_halo(fullOcean, mycells, 3)
@btime grow_halo(fullOcean, mycells, 3)

# @code_warntype divide_ocean(fullOcean, halowidth, 2, 2)


myCells = cellsInChunk[rank] # cellsedgesvertices[1]
myEdges = edgesInChunk[rank] # cellsedgesvertices[2]
myVertices = verticesInChunk[rank] # cellsedgesvertices[3]

mpasOcean = mpas_subset(fullOcean, myCells, myEdges, myVertices)



# mpasOcean = fullOcean

lYedge = maximum(fullOcean.yEdge) - minimum(fullOcean.yEdge)
lateralProfilePeriodic(y) = 1e-3*cos(y/lYedge * 4 * pi)
kelvinWaveExactNV, kelvinWaveExactSSH, kelvinWaveExactSolution!, boundaryCondition! = kelvinWaveGenerator(mpasOcean, lateralProfilePeriodic)






exacttime = 0.0

@btime calculate_normal_velocity_tendency!(mpasOcean)
@btime update_normal_velocity_by_tendency!(mpasOcean)

@btime update_normal_velocity_by_tendency_loop!(mpasOcean)


@btime boundaryCondition!(mpasOcean, exacttime)

@btime calculate_ssh_tendency!(mpasOcean)
@btime update_ssh_by_tendency!(mpasOcean)

#
#


@code_warntype calculate_normal_velocity_tendency!(mpasOcean)
@code_warntype update_normal_velocity_by_tendency!(mpasOcean)

@code_warntype boundaryCondition!(mpasOcean, exacttime)

@code_warntype calculate_ssh_tendency!(mpasOcean)
