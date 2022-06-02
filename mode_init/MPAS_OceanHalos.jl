import NCDatasets
import CUDA
import Base.Threads
# import KernelAbstractions

# include("Namelist.jl")
include("fixAngleEdge.jl")
include("MPAS_Ocean.jl")

mutable struct MPAS_OceanHalos

    # fullOcean::MPAS_Ocean

    haloChunks::Array{MPAS_Ocean}

    # cellsWithinHalo::AbstractArray
    # cellsOfHalo::AbstractArray

    # localCellsWithinHalo::AbstractArray
    # localCellsOfHalo::AbstractArray

    cellsToRank::Array
    cellsFromRank::Array

    edgesToRank::Array
    edgesFromRank::Array

    # load state from file
    function MPAS_OceanHalos(mpasOcean, haloWidth, nXChunks, nYChunks)
       mpasOceanHalos = new([], [], [], [], [])

        # divide up cells of mpasOcean into rectangular grid

        innerCells = collect([ [] for i in 1:nXChunks, j in 1:nYChunks])
        mpasOceanHalos.cellsFromRank = collect([ [] for i in 1:length(innerCells), j in 1:length(innerCells)])
        mpasOceanHalos.cellsToRank = collect([ [] for i in 1:length(innerCells), j in 1:length(innerCells)])

        chunkWidth = mpasOcean.lX/nXChunks
        chunkHeight = mpasOcean.lY/nYChunks

        for iCell in 1:mpasOcean.nCells
            append!(innerCells[Int(ceil(mpasOcean.xCell[iCell] / chunkWidth)),
                               Int(ceil(mpasOcean.yCell[iCell] / chunkHeight))], [iCell])
        end

        # make sub ocean objects

        for (i, chunk) in enumerate(innerCells)
            halo = grow_halo(mpasOcean, chunk, haloWidth)

            cells    = union(chunk, halo)
            edges    = Set(mpasOcean.edgesOnCell[cells])
            vertices = Set(mpasOcean.verticesOnCell[cells])

	    for (j, otherchunk) in enumerate(innerCells)
		mpasOceanHalos.cellsFromRank[i,j] = findall(iCell -> iCell in otherchunk, cells)
		mpasOceanHalos.cellsToRank[j,i]   = findall(iCell -> iCell in halo, otherchunk)
	    end


            append!(mpasOceanHalos.haloChunks, [mpas_subset(mpasOcean, cells, edges, vertices)])
        end

        return mpasOceanHalos
    end
end

function divide_ocean(mpasOcean::MPAS_Ocean, haloWidth, nXChunks, nYChunks)

    innerCells = collect([ [] for i in 1:nXChunks, j in 1:nYChunks])

	lXedge = maximum(mpasOcean.xEdge) - minimum(mpasOcean.xEdge)
	lYedge = maximum(mpasOcean.yEdge) - minimum(mpasOcean.yEdge)

	chunkWidth = lXedge/nXChunks
	chunkHeight = lYedge/nYChunks

    for iCell in 1:mpasOcean.nCells
        append!(innerCells[Int(ceil(mpasOcean.xCell[iCell] / chunkWidth)),
                           Int(ceil(mpasOcean.yCell[iCell] / chunkHeight))], [iCell])
    end

	cellsFromChunk = collect([ [] for i in 1:length(innerCells)])
	cellsToChunk = collect([ [] for i in 1:length(innerCells)])
	# cellsInChunk = []
	# edgesInChunk = []
	# verticesInChunk = []
	listtype = Vector{Any}
	haloCells = Array{listtype}(undef, length(innerCells))
	cellsInChunk = Array{listtype}(undef, length(innerCells))
	edgesInChunk = Array{listtype}(undef, length(innerCells))
	verticesInChunk = Array{listtype}(undef, length(innerCells))
	Threads.@threads for i in 1:length(innerCells)
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


	Threads.@threads for i in 1:length(innerCells) # (iChunk == "all" ? enumerate(innerCells) : [(1, innerCells[iChunk])])
	    chunk = innerCells[i]
	    cells = union(chunk, haloCells[i])
	    for (j, otherchunk) in enumerate(innerCells)
		    if j != i
				# if local cell (should be in the halo) overlaps the main (non-halo) part of another chunk, pull data from that chunk to update the halo
				localcells = findall(iCell -> iCell in otherchunk, cells)
				if length(localcells) > 0
					order = sortperm(cellsInChunk[i][localcells]) # sort relative to global index so local index order is consistent
					append!(cellsFromChunk[i], [(j, localcells[order])])
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

	return cellsInChunk, edgesInChunk, verticesInChunk, cellsFromChunk, cellsToChunk
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



function grow_halo(mpasOcean, cells, radius)
    # grow halo from edge of region
    outerCells = [cells...]
    newCells = []
    haloCells = []

    for i in 1:radius
        for iCell in outerCells
	    	if iCell != 0
                for jCell in mpasOcean.cellsOnCell[:,iCell]
					if ! (jCell in outerCells) && !(jCell in cells) && !(jCell in haloCells) && !(jCell == 0)
                        append!(newCells, jCell)
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


    function indexMap(iGlobal, indices)
		iLocals = findall(ind -> ind==iGlobal, indices)
		if length(iLocals) == 0
		    return 0
		end
		return iLocals[1]
    end

    ## take subsection of given cells, and map index fields to new local indices

    for field in cellCenteredFields
        setfield!(mpasSubOcean, field, getfield(mpasSubOcean, field)[collect(cells)])
    end
    for field in cell2dFields
		setfield!(mpasSubOcean, field, getfield(mpasSubOcean, field)[:,collect(cells)])
    end
    for field in cell2djiFields
		setfield!(mpasSubOcean, field, getfield(mpasSubOcean, field)[collect(cells),:])
    end
    for field in cellIndexFields
		setfield!(mpasSubOcean, field, map(i->indexMap(i,collect(cells)), getfield(mpasSubOcean, field)))
    end
    mpasSubOcean.nCells = length(cells)


    for field in edgeCenteredFields
        setfield!(mpasSubOcean, field, getfield(mpasSubOcean, field)[collect(edges)])
    end
    for field in edge2dFields
		setfield!(mpasSubOcean, field, getfield(mpasSubOcean, field)[:,collect(edges)])
    end
	for field in edge2djiFields
		setfield!(mpasSubOcean, field, getfield(mpasSubOcean, field)[collect(edges),:])
    end
    for field in edgeIndexFields
		setfield!(mpasSubOcean, field, map(i->indexMap(i,collect(edges)), getfield(mpasSubOcean, field)))
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
		setfield!(mpasSubOcean, field, map(i->indexMap(i,collect(vertices)), getfield(mpasSubOcean, field)))
    end
    mpasSubOcean.nVertices = length(vertices)

    return mpasSubOcean

end


function moveArrays!(mpasOceanHalos::MPAS_OceanHalos, array_type)
    for mpasOcean in mpasOceanHalos.HaloChunks
        moveArrays!(mpasOcean, array_type)
    end
end
