

using TimerOutputs
using BenchmarkTools
using DelimitedFiles

CODE_ROOT = pwd() * "/"#../"

include(CODE_ROOT * "mode_init/MPAS_Ocean.jl")
include(CODE_ROOT * "mode_init/MPAS_OceanHalos.jl")

include(CODE_ROOT * "mode_init/initial_conditions.jl")
include(CODE_ROOT * "mode_forward/time_steppers.jl")
include(CODE_ROOT * "mode_init/exactsolutions.jl")

rank = 1
worldrank = 1
root = 1
nprocs = 4
nCellsX = 64
halowidth = 5

function rootprint(text)
	if rank == root
		println(text)
	end
end

# println("loading ocean from file")
fullOcean =  MPAS_Ocean(
				CODE_ROOT * "MPAS_O_Shallow_Water/MPAS_O_Shallow_Water_Mesh_Generation/CoastalKelvinWaveMesh/ConvergenceStudyMeshes",
		       "culled_mesh_$nCellsX.nc", "mesh_$nCellsX.nc", nvlevels=1, periodicity="NonPeriodic_x")



partitiondir = CODE_ROOT * "scaling_test/graphparts/$(nCellsX)x$(nCellsX)/"
# fpath = CODE_ROOT * "output/kelvinwave/resolution$(nCellsX)x$(nCellsX)/steps$(halowidth*ncycles)"


println("dividing ocean")
partitionfile = partitiondir * "graph.info.part.$nprocs"
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

cellsOnCell = replace(readdlm(partitiondir * "graph.info"), ""=>0)

function partition(partitionfile::String = partitiondir * "graph.info.part.$nprocs")
	global mycells, haloCells
	# rootprint("partitioning $nCellsX x $nCellsX x $(fullOcean.nVertLevels) ocean among $nprocs processes with $(partitionfile)")
	partitions::Array{Int64,1} = dropdims((readdlm(partitionfile)), dims=2)
	cellsOnCell::Array{Int64,2} = replace(readdlm(partitiondir * "graph.info"), ""=>0)
	# partitions = Array{Int64,1}()
	# println("arr size ", sizeof(partitions))
	mycells = findall(proc -> proc == worldrank-1, partitions)




	# rootprint("halos: $halowidth wide")
	haloCells = grow_halo(cellsOnCell, mycells, halowidth)
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

# @code_warntype partition()

# @btime partition()

# @btime partitionorig()

partition()


# @code_warntype divide_ocean(fullOcean, halowidth, 2, 2)

function subset()
	myCells     = union(mycells, haloCells)
	mpasOcean =  MPAS_Ocean(
					CODE_ROOT * "MPAS_O_Shallow_Water/MPAS_O_Shallow_Water_Mesh_Generation/CoastalKelvinWaveMesh/ConvergenceStudyMeshes",
			       "culled_mesh_$nCellsX.nc", "mesh_$nCellsX.nc", nvlevels=1, periodicity="NonPeriodic_x",
				   cells=:)
end
# @code_warntype subset()
# @btime subset()
subset()


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
