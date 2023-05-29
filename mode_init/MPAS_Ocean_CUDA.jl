import CUDA
using SparseArrays

"""
    GPU-located copy of the prognostic fields and arrays used in computation to advance prog. fields.
    Does not contain everything MPAS_Ocean does, only these fields.
"""
mutable struct MPAS_Ocean_CUDA

    # cell centered fields
#     sshTendency::CUDA.CuArray{Float64,1}
    # sshCurrent::CUDA.CuArray{Float64,1}
    layerThicknessTendency::CUDA.CuArray{Float64,2}
    layerThickness::CUDA.CuArray{Float64,2}
    bottomDepth::CUDA.CuArray{Float64,1}
    areaCell::CUDA.CuArray{Float64,1}

    nEdgesOnCell::CUDA.CuArray{Int64,1}

    edgesOnCell::CUDA.CuArray{Int64,2}
    cellsOnCell::CUDA.CuArray{Int64,2}

    edgeSignOnCell::CUDA.CuArray{Int64,2}
    
    maxLevelCell::CUDA.CuArray{Int64,1}

    # edge centered fields
    fEdge::CUDA.CuArray{Float64,1}
    dcEdge::CUDA.CuArray{Float64,1}
    dvEdge::CUDA.CuArray{Float64,1}
    yEdge::CUDA.CuArray{Float64,1}
    xEdge::CUDA.CuArray{Float64,1}
    angleEdge::CUDA.CuArray{Float64,1}
    
    maxLevelEdgeTop::CUDA.CuArray{Int64,1}

    
    normalVelocityTendency::CUDA.CuArray{Float64,2}
    normalVelocityCurrent::CUDA.CuArray{Float64,2}
    weightsOnEdge::CUDA.CuArray{Float64,2}


    nEdgesOnEdge::CUDA.CuArray{Int64,1}

    cellsOnEdge::CUDA.CuArray{Int64,2}
    edgesOnEdge::CUDA.CuArray{Int64,2}

    boundaryEdge::CUDA.CuArray{Int64,2} # only 0 or 1, could be less exact datatype. But will that slow down more by causing casting?

    # scalars
    nCells::Int64
    nEdges::Int64

    gravity::Float64
    dt::Float64


    function MPAS_Ocean_CUDA(mpasOcean::MPAS_Ocean)# where Float64<:AbstractFloat

        mpasOceanCuda = new()

        for (field, typ) in zip(fieldnames(MPAS_Ocean_CUDA), MPAS_Ocean_CUDA.types)
            setfield!(mpasOceanCuda, field, convert(typ, getfield(mpasOcean, field)))
        end

        return mpasOceanCuda
    end
end
