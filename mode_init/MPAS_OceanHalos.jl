import NCDatasets
import CUDA
import Base.Threads
# import KernelAbstractions

# include("Namelist.jl")
include("fixAngleEdge.jl")
include("MPAS_Ocean.jl")

mutable struct MPAS_OceanHalos

    # list of indices in main MPAS ocean which this rank is managing
    myCells::Array{Int64}

    # maps <rank index +1> to a list of cell indices that rank has this rank needs, or vice versa
    cellsToMyRankFrom::Dict{Int64, Array}
    cellsFromMyRankTo::Dict{Int64, Array}

    # similar to cell Dicts, but with edges
    edgesToMyRankFrom::Dict{Int64, Array}
    edgesFromMyRankTo::Dict{Int64, Array}

    # maps <rank index +1> to an array equal in size and shape to the arrays of the prognostic variables of cells in the halo
    haloBuffersCells::Dict{Int64, Array}
    haloBuffersEdges::Dict{Int64, Array}
    

    # load state from file
    function MPAS_OceanHalos(comm, cellsOnCell, partitions)
        haloMesh = new()

	mycells = findall(proc -> proc == worldrank-1, partitions)

	(haloCells = grow_halo(cellsOnCell, mycells, halowidth))

	cellsFromRank = Dict{Int64, Array}()
	for iCell in haloCells
	    chunk = partitions[iCell] + 1
	    if chunk in keys(cellsFromRank)
		push!(cellsFromRank[chunk], iCell)
	    else
		cellsFromRank[chunk] = [iCell]
	    end
	end

	# tell neighbors which cells they should send me
	cellsToRank = Dict{Int64, Array}()
	tag = 2 # tag for communicating halo mesh information
	sendreqs = Array{MPI.Request,1}()
	recvreqs = Array{MPI.Request,1}()

	MPI.Barrier(comm)
	for (chunk, cells) in cellsFromRank
	    push!(sendreqs, MPI.Isend(cells, chunk-1, tag, comm))
	end
	MPI.Barrier(comm)
	for (chunk, cells) in cellsFromRank
	    numcells = MPI.Get_count(MPI.Probe(chunk-1, tag, comm), Int64)
	    cellsToRank[chunk] = Array{Int64,1}(undef, numcells)
	    push!(recvreqs, MPI.Irecv!(cellsToRank[chunk], chunk-1, tag, comm))
	end
	MPI.Waitall!(vcat(sendreqs, recvreqs))

	haloMesh.myCells     = union(mycells, haloCells)

	# order the arrays based on the local index they will have
	haloMesh.cellsToMyRankFrom = Dict{Int64, Array}()
	for (chunk, cells) in cellsToRank
	    localcells = globalToLocal(cells, haloMesh.myCells)
	    if length(localcells) > 0
		order = sortperm(haloMesh.myCells[localcells]) # order according to global index, so received in correct place
		haloMesh.cellsToMyRankFrom[chunk] = localcells[order]
	    end
	end
	haloMesh.cellsFromMyRankTo = Dict{Int64, Array}()
	for (chunk, cells) in cellsFromRank
	    localcells = globalToLocal(cells, haloMesh.myCells)
	    if length(localcells) > 0
		order = sortperm(haloMesh.myCells[localcells]) # order according to global index, so received in correct place
		haloMesh.cellsFromMyRankTo[chunk] = localcells[order]
	    end
	end



	return haloMesh

#
    end
end

function fillEdgeArraysFromCells!(haloMesh, edgesOnCell)
	myEdges = collect(Set(edgesOnCell[:,:])) # ocean will already be subsetted, so the cells are all the cells

	haloMesh.edgesToMyRankFrom = Dict{Int64, Array}()
	for (chunk, localcells) in haloMesh.cellsToMyRankFrom
	    localedges = collect(Set(edgesOnCell[:,localcells])) # Set() to remove duplicates
	    if length(localedges) > 0
		order = sortperm(myEdges[localedges])  # order according to global index, so received in correct place
		haloMesh.edgesToMyRankFrom[chunk] = localedges[order]
	    end
	end
	haloMesh.edgesFromMyRankTo = Dict{Int64, Array}()
	for (chunk, localcells) in haloMesh.cellsFromMyRankTo
	    localedges = collect(Set(edgesOnCell[:,localcells])) # Set() to remove duplicates
	    if length(localedges) > 0
		order = sortperm(myEdges[localedges])  # order according to global index, so received in correct place
		haloMesh.edgesFromMyRankTo[chunk] = localedges[order]
	    end
	end
end

function generateBufferArrays!(haloMesh, mpasOcean)

	haloMesh.haloBuffersCells = Dict{Int64, Array}() # temporarily stores new halo layer thickness (so all can be recieved before overwriting)
	for (srcrank, localcells) in (haloMesh.cellsToMyRankFrom)
		haloMesh.haloBuffersCells[srcrank] = Array{eltype(mpasOcean.layerThickness)}(undef, mpasOcean.nVertLevels, length(localcells))
	end
	haloMesh.haloBuffersEdges = Dict{Int64, Array}() # temporarily stores new halo normal velocity
	for (srcrank, localedges) in (haloMesh.edgesToMyRankFrom)
		haloMesh.haloBuffersEdges[srcrank] = Array{eltype(mpasOcean.normalVelocityCurrent)}(undef, mpasOcean.nVertLevels, length(localedges)) 
	end

end







function divide_ocean(mpasOcean::MPAS_Ocean, haloWidth, nXChunks, nYChunks)

#        # divide up cells of mpasOcean into rectangular grid
#
#        innerCells = collect([ [] for i in 1:nXChunks, j in 1:nYChunks])
#        haloMesh.cellsFromRank = collect([ [] for i in 1:length(innerCells), j in 1:length(innerCells)])
#        haloMesh.cellsToRank = collect([ [] for i in 1:length(innerCells), j in 1:length(innerCells)])
#
#        chunkWidth = mpasOcean.lX/nXChunks
#        chunkHeight = mpasOcean.lY/nYChunks
#
#        for iCell in 1:mpasOcean.nCells
#            append!(innerCells[Int(ceil(mpasOcean.xCell[iCell] / chunkWidth)),
#                               Int(ceil(mpasOcean.yCell[iCell] / chunkHeight))], [iCell])
#        end
#
#        # make sub ocean objects
#
#        for (i, chunk) in enumerate(innerCells)
#            halo = grow_halo(mpasOcean, chunk, haloWidth)
#
#            cells    = union(chunk, halo)
#            edges    = Set(mpasOcean.edgesOnCell[cells])
#            vertices = Set(mpasOcean.verticesOnCell[cells])
#
#	    for (j, otherchunk) in enumerate(innerCells)
#		haloMesh.cellsFromRank[i,j] = findall(iCell -> iCell in otherchunk, cells)
#		haloMesh.cellsToRank[j,i]   = findall(iCell -> iCell in halo, otherchunk)
#	    end
#
#
#            append!(haloMesh.haloChunks, [mpas_subset(mpasOcean, cells, edges, vertices)])
#        end
    innerCells = collect([ [] for i in 1:nXChunks, j in 1:nYChunks])

	lXedge = maximum(mpasOcean.xEdge) - minimum(mpasOcean.xEdge)
	lYedge = maximum(mpasOcean.yEdge) - minimum(mpasOcean.yEdge)

	chunkWidth = lXedge/nXChunks
	chunkHeight = lYedge/nYChunks

    for iCell in 1:mpasOcean.nCells
        append!(innerCells[Int(ceil(mpasOcean.xCell[iCell] / chunkWidth)),
                           Int(ceil(mpasOcean.yCell[iCell] / chunkHeight))], [iCell])
    end

	haloMesh.cellsFromRank = collect([ [] for i in 1:length(innerCells)])
	cellsToChunk = collect([ [] for i in 1:length(innerCells)])
	# cellsInChunk = []
	# edgesInChunk = []
	# verticesInChunk = []
	listtype = Vector{Any}
	haloCells = Array{listtype}(undef, length(innerCells))
	cellsInChunk = Array{listtype}(undef, length(innerCells))
	edgesInChunk = Array{listtype}(undef, length(innerCells))
	verticesInChunk = Array{listtype}(undef, length(innerCells))
	# Threads.@threads
	for i in 1:length(innerCells)
	    chunk = innerCells[i]
        halo = grow_halo(mpasOcean, chunk, haloWidth)
	    haloCells[i] = halo

        cells    = union(chunk, halo)
        edges    = collect(Set(mpasOcean.edgesOnCell[:,cells]))
        vertices = collect(Set(mpasOcean.verticesOnCell[:,cells]))

	    cellsInChunk[i] = cells
	    edgesInChunk[i] = edges
	    verticesInChunk[i] = vertices

    end


	# Threads.@threads
	for i in 1:length(innerCells) # (iChunk == "all" ? enumerate(innerCells) : [(1, innerCells[iChunk])])
	    chunk = innerCells[i]
	    cells = union(chunk, haloCells[i])
	    for (j, otherchunk) in enumerate(innerCells)
		    if j != i
				# if local cell (should be in the halo) overlaps the main (non-halo) part of another chunk, pull data from that chunk to update the halo
				localcells = findall(iCell -> iCell in otherchunk, cells)
				if length(localcells) > 0
					order = sortperm(cellsInChunk[i][localcells]) # sort relative to global index so local index order is consistent
					append!(haloMesh.cellsFromRank[i], [(j, localcells[order])])
				end
		    end
	    end
	    for (j, otherhalo) in enumerate(haloCells)
		    if j != i
				# if local cell in main (non-halo) area overlaps another chunk's halo, send local data to that chunk to update its halo
				localcells = findall(iCell -> iCell in otherhalo, chunk)
				if length(localcells) > 0
					order = sortperm(cellsInChunk[i][localcells])
					append!(cellsToChunk[i], [(j, localcells[order])])
				end
		    end
	    end
	end

	return haloMesh #cellsInChunk, edgesInChunk, verticesInChunk, haloMesh.cellsFromRank, cellsToChunk
end


#         function make_region(startCell)
#             # grow region out from start cell

#             innerCells = grow_halo(mpasOcean, [startCell], 20)

#             append!(innerCells, outerCells)

#             haloCells = grow_halo(mpasOcean, innerCells, halo_width)

#             append!(haloCells, outerCells)

#             return innerCells, haloCells
#         end


sizeof_main_arrays(mpasOcean::MPAS_Ocean) = sizeof_main_arrays(mpasOcean, collect(1:mpasOcean.nCells), collect(1:mpasOcean.nEdges))
function sizeof_main_arrays(mpasOcean::MPAS_Ocean, cells, edges)
    mainCellArrayNames = [
        :sshTendency,
        :sshCurrent,
        :nEdgesOnCell,
        :edgesOnCell,
        :cellsOnCell,
        :bottomDepth,
        :areaCell,
        :edgeSignOnCell,
    ]

    mainEdgeArrayNames = [
        :normalVelocityTendency,
        :normalVelocityCurrent,
        :cellsOnEdge,
        :nEdgesOnEdge,
        :edgesOnEdge,
        :weightsOnEdge,
        :fEdge,
        :dcEdge,
        :dvEdge,
    ]

    size::Integer = 0

    for name in mainCellArrayNames
        size += sizeof(getfield(mpasOcean, name)[cells])
    end
    for name in mainEdgeArrayNames
        size += sizeof(getfield(mpasOcean, name)[edges])
    end

    return size
end

function grow_halo(cellsOnCell::Array{Int64,2}, cells::Array{Int64,1}, radius::Int64)
    # grow halo from edge of region
    outerCells::Array{Int64,1} = copy(cells)
    newCells::Array{Int64,1} = []
    haloCells::Array{Int64,1} = []

    for i in 1:radius
        for iCell in outerCells
	    	if iCell != 0
                for jCell in cellsOnCell[iCell,:]
					if !(jCell in newCells) && !(jCell in outerCells) && !(jCell in cells) && !(jCell in haloCells) && !(jCell == 0)
                        append!(newCells, jCell) #
                    end
				end
            end
        end
		append!(haloCells, newCells)
        outerCells = newCells
        newCells = []
    end

    return haloCells
end

function grow_halo(mpasOcean::MPAS_Ocean, cells::Array{Int64,1}, radius::Int64)
    return grow_halo(mpasOcean.cellsOnCell, cells, radius)
end

function mpas_subset(mpasOcean::MPAS_Ocean, cells, edges, vertices)
    mpasSubOcean = deepcopy(mpasOcean)

    cellCenteredFields = [
        :sshCurrent,
        :sshTendency,
        :bottomDepth,
        :latCell,
        :lonCell,
        :xCell,
        :yCell,
        :areaCell,
        :fCell,
        :maxLevelCell,
        :gridSpacing,
        :nEdgesOnCell,
    ]
    cell2dFields = [
        :cellsOnCell,
        :edgesOnCell,
        :verticesOnCell,
    ]
    cell2djiFields = [
		:edgeSignOnCell,
        :boundaryCell,
		:cellMask,
    ]
    cellIndexFields = [
		:cellsOnCell,
		:cellsOnEdge,
		:cellsOnVertex
    ]



    edgeCenteredFields = [
        :xEdge,
        :yEdge,
        :dvEdge,
        :dcEdge,
        :fEdge,
        :angleEdge,
        :maxLevelEdgeTop,
        :maxLevelEdgeBot,
        :nEdgesOnEdge,
    ]
    edge2dFields = [
        :cellsOnEdge,
        :edgesOnEdge,
        :verticesOnEdge,
		:weightsOnEdge,
    ]
    edge2djiFields = [
        :normalVelocityCurrent,
        :normalVelocityTendency,
		:boundaryEdge,
		:edgeMask,
    ]
    edgeIndexFields = [
		:edgesOnCell,
		:edgesOnEdge,
		:edgesOnVertex,
    ]




    vertexCenteredFields = [
        :latVertex,
        :lonVertex,
        :xVertex,
        :yVertex,
        :fVertex,
        :areaTriangle,
        :maxLevelVertexTop,
        :maxLevelVertexBot,
    ]
    vertex2dFields = [
        :cellsOnVertex,
        :edgesOnVertex,
    ]
    vertex2djiFields = [
		:edgeSignOnVertex,
        :boundaryVertex,
		:vertexMask
    ]
    vertexIndexFields = [
		:verticesOnCell,
		:verticesOnEdge,
    ]


    mycells = collect(cells)
    globtolocalcell = Dict{Int64,Int64}()
    for (iLocal, iGlobal) in enumerate(mycells)
	    globtolocalcell[iGlobal] = iLocal
    end
    myedges = collect(edges)
    globtolocaledge = Dict{Int64,Int64}()
    for (iLocal, iGlobal) in enumerate(myedges)
	    globtolocaledge[iGlobal] = iLocal
    end
    myvertices = collect(vertices)
    globtolocalvertex = Dict{Int64,Int64}()
    for (iLocal, iGlobal) in enumerate(myvertices)
	    globtolocalvertex[iGlobal] = iLocal
    end
	#=
    function indexMap(iGlobal, indices)
		iLocals = findall(ind -> ind==iGlobal, indices)
		if length(iLocals) == 0
		    return 0
		end
		return iLocals[1]
    end
    =#
    #=
	function globalToLocal(globalcells, mycells)
		localcells = findall(iCell -> iCell in globalcells, mycells)
		return localcells
	end
	=#
	function globalToLocal(globals, globmap)
		map(iGlobal -> (if iGlobal in keys(globmap) return globmap[iGlobal] else return 0; end),
		    	globals)
	end

    ## take subsection of given cells, and map index fields to new local indices

    for field in cellCenteredFields
        setfield!(mpasSubOcean, field, getfield(mpasSubOcean, field)[mycells])
    end
    for field in cell2dFields
		setfield!(mpasSubOcean, field, getfield(mpasSubOcean, field)[:,mycells])
    end
    for field in cell2djiFields
		setfield!(mpasSubOcean, field, getfield(mpasSubOcean, field)[mycells,:])
    end
    for field in cellIndexFields
		setfield!(mpasSubOcean, field, # map(i->indexMap(i,collect(cells)),
			  globalToLocal(getfield(mpasSubOcean, field), globtolocalcell))
    end
    mpasSubOcean.nCells = length(cells)

    for field in edgeCenteredFields
        setfield!(mpasSubOcean, field, getfield(mpasSubOcean, field)[myedges])
    end
    for field in edge2dFields
		setfield!(mpasSubOcean, field, getfield(mpasSubOcean, field)[:,myedges])
    end
	for field in edge2djiFields
		setfield!(mpasSubOcean, field, getfield(mpasSubOcean, field)[myedges,:])
    end
    for field in edgeIndexFields
	    setfield!(mpasSubOcean, field, globalToLocal(getfield(mpasSubOcean, field), globtolocaledge))
    end
    mpasSubOcean.nEdges = length(edges)




    for field in vertexCenteredFields
        setfield!(mpasSubOcean, field, getfield(mpasSubOcean, field)[collect(vertices)])
    end
    for field in vertex2dFields
		setfield!(mpasSubOcean, field, getfield(mpasSubOcean, field)[:,collect(vertices)])
    end
    for field in vertex2djiFields
		setfield!(mpasSubOcean, field, getfield(mpasSubOcean, field)[collect(vertices),:])
    end
    for field in vertexIndexFields
		setfield!(mpasSubOcean, field, globalToLocal(getfield(mpasSubOcean, field), globtolocalvertex))
    end
    mpasSubOcean.nVertices = length(vertices)


    return mpasSubOcean

end
