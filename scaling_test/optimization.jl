using MPI
MPI.Init()

comm = MPI.COMM_WORLD
worldrank = MPI.Comm_rank(comm) + 1
root = 1
commsize = MPI.Comm_size(comm)

using TimerOutputs
using BenchmarkTools
using DelimitedFiles

CODE_ROOT = pwd() * "/"#../"

include(CODE_ROOT * "mode_init/MPAS_Ocean.jl")
include(CODE_ROOT * "mode_init/MPAS_OceanHalos.jl")

include(CODE_ROOT * "mode_init/exactsolutions.jl")
include(CODE_ROOT * "mode_forward/time_steppers.jl")

include(CODE_ROOT * "mode_forward/update_halos.jl")


rank = 1
worldrank = 1
root = 1
nprocs = 4
nCellsX = 128
halowidth = 1

nvlevels = 100

function rootprint(text)
	if rank == root
		println(text)
	end
end

meshpath = CODE_ROOT * "MPAS_O_Shallow_Water/MPAS_O_Shallow_Water_Mesh_Generation/CoastalKelvinWaveMesh/ConvergenceStudyMeshes"


partitiondir = CODE_ROOT * "scaling_test/graphparts/$(nCellsX)x$(nCellsX)/"
# fpath = CODE_ROOT * "output/kelvinwave/resolution$(nCellsX)x$(nCellsX)/steps$(halowidth*ncycles)"


function globalToLocal(globalcells, mycells)
	localcells = findall(iCell -> iCell in globalcells, mycells)
	return localcells
end



# @code_warntype divide_ocean(fullOcean, halowidth, 2, 2)


# mpasOcean = fullOcean

# lYedge = maximum(fullOcean.yEdge) - minimum(fullOcean.yEdge)
# lateralProfilePeriodic(y) = 1e-3*cos(y/lYedge * 4 * pi)
# kelvinWaveExactNV, kelvinWaveExactSSH, kelvinWaveExactSolution!, boundaryCondition! = kelvinWaveGenerator(mpasOcean, lateralProfilePeriodic)





exacttime = 0.0

# @btime calculate_normal_velocity_tendency!(mpasOcean)
# @btime update_normal_velocity_by_tendency!(mpasOcean)

# @btime update_normal_velocity_by_tendency_loop!(mpasOcean)

# @btime boundaryCondition!(mpasOcean, exacttime)

# @btime calculate_ssh_tendency!(mpasOcean)
# @btime update_ssh_by_tendency!(mpasOcean)


#
#


# @code_warntype calculate_normal_velocity_tendency!(mpasOcean)
# @code_warntype update_normal_velocity_by_tendency!(mpasOcean)

# @code_warntype boundaryCondition!(mpasOcean, exacttime)

# @code_warntype calculate_ssh_tendency!(mpasOcean)


worldrank = rank
nprocs = 4

splitcomm = MPI.Comm_split(comm, Int(worldrank>nprocs), worldrank)






	meshpath = CODE_ROOT * "MPAS_O_Shallow_Water/MPAS_O_Shallow_Water_Mesh_Generation/CoastalKelvinWaveMesh/ConvergenceStudyMeshes"

	# kelvin wave initial condition and exact solution functions
	if worldrank == root
		fullOcean = MPAS_Ocean(meshpath, "culled_mesh_$nCellsX.nc", "mesh_$nCellsX.nc", periodicity="NonPeriodic_x", nvlevels=nvlevels)
		lYedge = maximum(fullOcean.yEdge) - minimum(fullOcean.yEdge)
		meanCoriolisParameterf = sum(fullOcean.fEdge) / length(fullOcean.fEdge)
		meanFluidThicknessH = sum(fullOcean.bottomDepth)/length(fullOcean.bottomDepth)
		c = sqrt(fullOcean.gravity*meanFluidThicknessH)
		rossbyRadiusR = c/meanCoriolisParameterf
		gravity = fullOcean.gravity
		dt = fullOcean.dt
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


partitionfile = partitiondir * "graph.info.part.$nprocs"
rootprint("partitioning $nCellsX x $nCellsX x $(nvlevels) ocean among $nprocs processes with $(partitionfile)")
partitions = dropdims(readdlm(partitionfile), dims=2)
mycells = findall(proc -> proc == worldrank-1, partitions)

	cellsOnCell = Array{Int64,2}(replace(readdlm(partitiondir * "graph.info")[2:end,:], ""=>0))
(haloCells = grow_halo(cellsOnCell, mycells, halowidth))

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
sendreqs = Array{MPI.Request,1}()
recvreqs = Array{MPI.Request,1}()

MPI.Barrier(splitcomm)
for (chunk, cells) in cellsFromChunk
    push!(sendreqs, MPI.Isend(cells, chunk-1, tag, splitcomm))
end
MPI.Barrier(splitcomm)
for (chunk, cells) in cellsFromChunk
    numcells = MPI.Get_count(MPI.Probe(chunk-1, tag, splitcomm), Int64)
    cellsToChunk[chunk] = Array{Int64,1}(undef, numcells)
    push!(recvreqs, MPI.Irecv!(cellsToChunk[chunk], chunk-1, tag, splitcomm))
end
MPI.Waitall!(vcat(sendreqs, recvreqs))

rootprint("subsetting ocean")

#			@timeit timout "subset" begin
    myCells     = union(mycells, haloCells)
    rootprint("cell count $(length(myCells))")
    mpasOcean = MPAS_Ocean(meshpath, "culled_mesh_$nCellsX.nc", "mesh_$nCellsX.nc", periodicity="NonPeriodic_x", nvlevels=nvlevels, cells=myCells)
    mpasOcean.dt = dt #1e-2 * minimum(mpasOcean.dcEdge) / c
    rootprint("dt: $(mpasOcean.dt)")

    # rootprint("edgesoncell $(mpasOcean.edgesOnCell[1:10]")
#			end
    myEdges = collect(Set(mpasOcean.edgesOnCell[:,:]))

cellsToMyChunk = Dict{Int64, Array}()
for (chunk, cells) in cellsToChunk
    localcells = globalToLocal(cells, myCells)
    if length(localcells) > 0
        order = sortperm(myCells[localcells]) # order according to global index, so received in correct place
        cellsToMyChunk[chunk] = localcells[order]
    end
end
cellsFromMyChunk = Dict{Int64, Array}()
for (chunk, cells) in cellsFromChunk
    localcells = globalToLocal(cells, myCells)
    if length(localcells) > 0
        order = sortperm(myCells[localcells]) # order according to global index, so received in correct place
        cellsFromMyChunk[chunk] = localcells[order]
    end
end
edgesToMyChunk = Dict{Int64, Array}()
for (chunk, localcells) in cellsToMyChunk
    localedges = collect(Set(mpasOcean.edgesOnCell[:,localcells])) # Set() to remove duplicates
    if length(localedges) > 0
        order = sortperm(myEdges[localedges])  # order according to global index, so received in correct place
        edgesToMyChunk[chunk] = localedges[order]
    end
end
edgesFromMyChunk = Dict{Int64, Array}()
for (chunk, localcells) in cellsFromMyChunk
    localedges = collect(Set(mpasOcean.edgesOnCell[:,localcells])) # Set() to remove duplicates
    if length(localedges) > 0
        order = sortperm(myEdges[localedges])  # order according to global index, so received in correct place
        edgesFromMyChunk[chunk] = localedges[order]
    end
end





if rank == root
    @btime update_halos!(splitcomm, mpasOcean, cellsFromMyChunk, cellsToMyChunk, edgesFromMyChunk, edgesToMyChunk)
else
    update_halos!(splitcomm, mpasOcean, cellsFromMyChunk, cellsToMyChunk, edgesFromMyChunk, edgesToMyChunk)
end


if rank == root
    @code_warntype update_halos!(splitcomm, mpasOcean, cellsFromMyChunk, cellsToMyChunk, edgesFromMyChunk, edgesToMyChunk)
else
    update_halos!(splitcomm, mpasOcean, cellsFromMyChunk, cellsToMyChunk, edgesFromMyChunk, edgesToMyChunk)
end

