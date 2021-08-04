import NCDatasets
import CUDA
# import KernelAbstractions

# include("Namelist.jl")
include("fixAngleEdge.jl")
include("BoundaryCondition.jl")
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
        
        chunkCells = collect([ [] for i in 1:nXChunks, j in 1:nYChunks])
        mpasOceanHalos.cellsFromRank = collect([ [] for i in 1:length(chunkCells), j in 1:length(chunkCells)])
        mpasOceanHalos.cellsToRank = collect([ [] for i in 1:length(chunkCells), j in 1:length(chunkCells)])
        
        chunkWidth = mpasOcean.lX/nXChunks
        chunkHeight = mpasOcean.lY/nYChunks
        
        for iCell in 1:mpasOcean.nCells
            append!(chunkCells[Int(ceil(mpasOcean.xCell[iCell] / chunkWidth)),
                               Int(ceil(mpasOcean.yCell[iCell] / chunkHeight))], [iCell])
        end
        
        # make sub ocean objects
        
        for (i, chunk) in enumerate(chunkCells)
            halo = grow_halo(mpasOcean, chunk, haloWidth)

            cells    = union(chunk, halo)
            edges    = Set(mpasOcean.edgesOnCell[cells])
            vertices = Set(mpasOcean.verticesOnCell[cells])

	    for (j, otherchunk) in enumerate(chunkCells)
		mpasOceanHalos.cellsFromRank[i,j] = findall(iCell -> iCell in otherchunk, cells)
		mpasOceanHalos.cellsToRank[j,i]   = findall(iCell -> iCell in halo, otherchunk)
	    end

            
            append!(mpasOceanHalos.haloChunks, [mpas_subset(mpasOcean, cells, edges, vertices)])
        end
        
        return mpasOceanHalos
    end
end

function divide_ocean(mpasOcean::MPAS_Ocean, haloWidth, nXChunks, nYChunks)
	cellsInChunk = []
	edgesInChunk = []
	verticesInChunk = []

        chunkCells = collect([ [] for i in 1:nXChunks, j in 1:nYChunks])
	# if iChunk == "all"
		cellsFromChunk = collect([ [] for i in 1:length(chunkCells)])
		cellsToChunk = collect([ [] for i in 1:length(chunkCells)])
	# else
		# cellsFromChunk = [[]]
		# cellsToChunk = [[]]
	# end
        
        chunkWidth = mpasOcean.lX/nXChunks
        chunkHeight = mpasOcean.lY/nYChunks
        
        for iCell in 1:mpasOcean.nCells
            append!(chunkCells[Int(ceil(mpasOcean.xCell[iCell] / chunkWidth)),
                               Int(ceil(mpasOcean.yCell[iCell] / chunkHeight))], [iCell])
        end
        
        # put halo borders around chunks
	# if iChunk == "all"
		chunkHalos = collect([ [] for i in 1:length(chunkCells)]) 
	# else
		# chunkHalos = [[]]
	# end
        for (i, chunk) in enumerate(chunkCells) # (iChunk == "all" ? enumerate(chunkCells) : [(1, chunkCells[iChunk])])
            halo = grow_halo(mpasOcean, chunk, haloWidth)
	    chunkHalos[i] = halo

            cells    = union(chunk, halo)
            edges    = Set(mpasOcean.edgesOnCell[cells])
            vertices = Set(mpasOcean.verticesOnCell[cells])


	    append!(cellsInChunk, [cells])
	    append!(edgesInChunk, [edges])
	    append!(verticesInChunk, [vertices])
	    
        end
         
        for (i, chunk) in enumerate(chunkCells) # (iChunk == "all" ? enumerate(chunkCells) : [(1, chunkCells[iChunk])])
	    for (j, otherchunk) in enumerate(chunkCells)
	    	# if local cell (should be in the halo) overlaps the main (non-halo) part of another chunk, pull data from that chunk to update the halo
		append!(cellsFromChunk[i], (j, findall(iCell -> iCell in otherchunk, cells)))
	    for (j, otherhalo) in enumerate(chunkHalos)
	    	# if local cell in main (non-halo) area overlaps another chunk's halo, send local data to that chunk to update its halo
		append!(cellsToChunk[i],   (j, findall(iCell -> iCell in otherhalo, chunk)))
	    end
	end

	# if iChunk != "all"
		# cellsFromChunk = dropdims(cellsFromChunk, dims=1)
	# end

	return cellsInChunk, edgesInChunk, verticesInChunk, cellsFromChunk, cellsToChunk
end


#         function make_region(startCell)
#             # grow region out from start cell
            
#             chunkCells = grow_halo(mpasOcean, [startCell], 20)
            
#             append!(chunkCells, outerCells)
            
#             haloCells = grow_halo(mpasOcean, chunkCells, halo_width)
            
#             append!(haloCells, outerCells)
            
#             return chunkCells, haloCells
#         end


sizeof_main_arrays(mpasOcean::MPAS_Ocean) = sizeof_main_arrays(mpasOcean, collect(1:mpasOcean.nCells), collect(1:mpasOcean.nEdges))
function sizeof_main_arrays(mpasOcean::MPAS_Ocean, cells, edges)
    mainCellArrayNames = [
        :sshTendency,
        :ssh,
        :nEdgesOnCell,
        :edgesOnCell,
        :cellsOnCell,
        :bottomDepth,
        :areaCell,
        :edgeSignOnCell,
    ]
    
    mainEdgeArrayNames = [
        :normalVelocityTendency,
        :normalVelocity,
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
            for jCell in mpasOcean.cellsOnCell[iCell]
                if ! (jCell in outerCells) && !(jCell in cells) && !(jCell in haloCells)
                    append!(newCells, jCell)
		end
            end
        end
        if i > 1
		append!(haloCells, outerCells)
	end
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
        :edgeSignOnCell,
        :latCell,
        :lonCell,
        :xCell,
        :yCell,
        :areaCell,
        :fCell,
        :maxLevelCell,
        :gridSpacing,
        :boundaryCell,
    ]
    cellCenteredIndexFields = [
        :cellsOnCell,
        :edgesOnCell,
        :verticesOnCell,
        :kiteIndexOnCell,
        :nEdgesOnCell,
    ]

    edgeCenteredFields = [
        :normalVelocityCurrent,
        :normalVelocityTendency,
        # :nEdges,
        :xEdge,
        :yEdge,
        :dvEdge,
        :dcEdge,
        :fEdge,
        :angleEdge,
        :weightsOnEdge,
        :maxLevelEdgeTop,
        :maxLevelEdgeBot,
        :boundaryEdge,
        :edgeMask,
    ]
    edgeCenteredIndexFields = [
        :cellsOnEdge,
        :edgesOnEdge,
        :verticesOnEdge,
        :nEdgesOnEdge,
    ]
    
    vertexCenteredFields = [
        # :nVertices,
        :latVertex,
        :lonVertex,
        :xVertex,
        :yVertex,
        # :vertexDegree,
        :edgeSignOnVertex,
        :fVertex,
        :areaTriangle,
        :kiteAreasOnVertex,
        :maxLevelVertexTop,
        :maxLevelVertexBot,
        :boundaryVertex,
    ]
    vertexCenteredIndexFields = [
        :cellsOnVertex,
        :edgesOnVertex,
    ]

    function indexMap(iGlobal, indices)
	iLocals = findall(ind -> ind==iGlobal, indices)
	if length(iLocals) == 0
	    return 0
	end
	return iLocals[1]
    end

    for field in cellCenteredFields
        setfield!(mpasSubOcean, field, getfield(mpasSubOcean, field)[collect(cells)])
    end
    for field in cellCenteredIndexFields
	setfield!(mpasSubOcean, field, map(i -> indexMap(i, collect(cells)), getfield(mpasSubOcean, field)))
    end
    mpasSubOcean.nCells = length(cells)
    
    for field in edgeCenteredFields
        setfield!(mpasSubOcean, field, getfield(mpasSubOcean, field)[collect(edges)])
    end
    for field in edgeCenteredIndexFields
	setfield!(mpasSubOcean, field, map(i -> indexMap(i, collect(edges)), getfield(mpasSubOcean, field)))
    end
    mpasSubOcean.nEdges = length(edges)

    for field in vertexCenteredFields
        setfield!(mpasSubOcean, field, getfield(mpasSubOcean, field)[collect(vertices)])
    end
    for field in vertexCenteredIndexFields
	setfield!(mpasSubOcean, field, map(i -> indexMap(i, collect(vertices)), getfield(mpasSubOcean, field)))
    end
    mpasSubOcean.nVertices = length(vertices)

    return mpasSubOcean

end


function moveArrays!(mpasOceanHalos::MPAS_OceanHalos, array_type)
    for mpasOcean in mpasOceanHalos.HaloChunks
        moveArrays!(mpasOcean, array_type)
    end
end



