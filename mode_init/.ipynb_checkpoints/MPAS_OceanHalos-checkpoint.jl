import NCDatasets
import CUDA
# import KernelAbstractions

# include("Namelist.jl")
include("fixAngleEdge.jl")
include("BoundaryCondition.jl")


mutable struct MPAS_OceanHalos
    
    fullOcean::MPAS_Ocean
    
    haloChunks::Array{MPAS_Ocean}
    
    cellsWithinHalo::AbstractArray
    cellsOfHalo::AbstractArray
    
    
    # load state from file
    function MPAS_OceanHalos(mpasOcean, haloWidth, nXChunks, nYChunks)
       mpasOceanHalos = new(mpasOcean, [], [], [])
        
        # divide up cells of mpasOcean into rectangular grid
        
        chunkCells = collect([ [] for i in 1:nXChunks, j in 1:nYChunks])
        
        chunkWidth = mpasOcean.lX/nXChunks
        chunkHeight = mpasOcean.lY/nYChunks
        
        for iCell in 1:mpasOcean.nCells
            append!(chunkCells[Int(floor(mpasOcean.xCell[iCell] / chunkWidth)),
                               Int(floor(mpasOcean.yCell[iCell] / chunkHeight))], [iCell])
        end
        
        # make sub ocean objects
        
        for chunk in chunkCells
            append!(mpasOceanHalos.cellsWithinHalo, [chunk])
            halo = grow_halo(mpasOcean, chunk, haloWidth)
            append!(mpasOceanHalos.cellsOfHalo, [halo])
            
            cells    = union(chunk, halo)
            edges    = Set(mpasOceanHalos.fullOcean.edgesOnCell[cells])
            vertices = Set(mpasOceanHalos.fullOcean.verticesOnCell[cells])
            
            append!(mpasOceanHalos.haloChunks, mpas_subset(mpasOcean, cells, edges, vertices))
        end
        
        return mpasOceanHalos
    end
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

    haloCells = []
    outerCells = []
    newCells = []

    for i in 1:radius
        for iCell in outerCells
            for jCell in mpasOcean.cellsOnCell[iCell]
                if !jCell in outerCells && !jCell in cells && !jCell in haloCells
                    append!(newCells, jCell)
                end
            end
        end
        append!(haloCells, outerCells)
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
        :nCells,
        :cellsOnCell,
        :edgesOnCell,
        :verticesOnCell,
        :kiteIndexOnCell,
        :nEdgesOnCell,
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

    edgeCenteredFields = [
        :normalVelocityCurrent,
        :normalVelocityTendency,
        :nEdges,
        :cellsOnEdge,
        :edgesOnEdge,
        :verticesOnEdge,
        :nEdgesOnEdge,
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

    vertexCenteredFields = [
        :nVertices,
        :latVertex,
        :lonVertex,
        :xVertex,
        :yVertex,
        :vertexDegree,
        :cellsOnVertex,
        :edgesOnVertex,
        :edgeSignOnVertex,
        :fVertex,
        :areaTriangle,
        :kiteAreasOnVertex,
        :maxLevelVertexTop,
        :maxLevelVertexBot,
        :boundaryVertex,
    ]

    for field in cellCenteredFields
        setfield!(mpasSubOcean, field, getfield(mpasSubOcean, field)[cells])
    end
    for field in edgeCenteredFields
        setfield!(mpasSubOcean, field, getfield(mpasSubOcean, field)[edges])
    end
    for field in vertexCenteredFields
        setfield!(mpasSubOcean, field, getfield(mpasSubOcean, field)[vertices])
    end

    return mpasSubOcean

end


function moveArrays!(mpasOceanHalos::MPAS_OceanHalos, array_type)
    for mpasOcean in mpasOceanHalos.HaloChunks
        moveArrays!(mpasOcean, array_type)
    end
end



