include("../mode_init/MPAS_Ocean_CUDA.jl")

######## CUDA GPU tendency calculation methods

function calculate_normal_velocity_tendency_cuda!(mpasOcean::MPAS_Ocean_CUDA)
    CUDA.@cuda blocks=cld(mpasOcean.nEdges, 1024) threads=1024 maxregs=64 calculate_normal_velocity_tendency_cuda_kernel!(
                                                                        mpasOcean.nEdges,
                                                                        mpasOcean.normalVelocityTendency,
                                                                        mpasOcean.normalVelocityCurrent,
                                                                        mpasOcean.layerThickness,
                                                                        mpasOcean.maxLevelEdgeTop,
                                                                        mpasOcean.cellsOnEdge,
                                                                        mpasOcean.nEdgesOnEdge,
                                                                        mpasOcean.edgesOnEdge,
                                                                        mpasOcean.weightsOnEdge,
                                                                        mpasOcean.fEdge,
                                                                        mpasOcean.dcEdge,
                                                                        mpasOcean.gravity)
end

function calculate_normal_velocity_tendency_cuda_kernel!(nEdges,
                                                          normalVelocityTendency,
                                                          normalVelocity,
                                                          layerThickness,
                                                          maxLevelEdgeTop,
                                                          cellsOnEdge,
                                                          nEdgesOnEdge,
                                                          edgesOnEdge,
                                                          weightsOnEdge,
                                                          fEdge,
                                                          dcEdge,
                                                          gravity)

    iEdge = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
    if iEdge <= nEdges

        cell1 = cellsOnEdge[1,iEdge]
        cell2 = cellsOnEdge[2,iEdge]
        
        # gravity term: take gradient across edge
        for k in 1:maxLevelEdgeTop[iEdge]
            normalVelocityTendency[k,iEdge] = gravity * ( layerThickness[k,cell1] - layerThickness[k,cell2] ) / dcEdge[iEdge]
        end

        # coriolis term: weighted sum of norm. vel. of surrounding edges to approximate tangential vel.
        for i = 1:nEdgesOnEdge[iEdge]
            eoe = edgesOnEdge[i, iEdge]

            if eoe != 0
                for k in 1:maxLevelEdgeTop[iEdge]
                    normalVelocityTendency[k,iEdge] += weightsOnEdge[i,iEdge] * normalVelocity[k,eoe] * fEdge[eoe]
                end
            end
        end


    end
    return
end







function update_normal_velocity_by_tendency_cuda!(mpasOcean::MPAS_Ocean_CUDA)
    CUDA.@cuda blocks=cld(mpasOcean.nEdges, 1024) threads=1024 maxregs=64 update_normal_velocity_by_tendency_cuda_kernel!(
                                                                        mpasOcean.nEdges,
                                                                        mpasOcean.normalVelocityCurrent,
                                                                        mpasOcean.dt,
                                                                        mpasOcean.normalVelocityTendency,
                                                                        mpasOcean.maxLevelEdgeTop)
end

function update_normal_velocity_by_tendency_cuda_kernel!(nEdges,
                                                         normalVelocityCurrent,
                                                         dt,
                                                         normalVelocityTendency,
                                                         maxLevelEdgeTop)
    iEdge = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
    if iEdge <= nEdges
        
        for k in 1:maxLevelEdgeTop[iEdge]
            normalVelocityCurrent[k,iEdge] += dt * normalVelocityTendency[k,iEdge]
        end
        
    end
    return
end








function calculate_thickness_tendency_cuda!(mpasOcean::MPAS_Ocean_CUDA)
    CUDA.@cuda blocks=cld(mpasOcean.nCells, 1024) threads=1024 maxregs=64 calculate_thickness_tendency_cuda_kernel!(
                                                                        mpasOcean.nCells,
                                                                        mpasOcean.layerThicknessTendency,
                                                                        mpasOcean.layerThickness,
                                                                        mpasOcean.normalVelocityCurrent,
                                                                        mpasOcean.maxLevelEdgeTop,
                                                                        mpasOcean.nEdgesOnCell,
                                                                        mpasOcean.edgesOnCell,
                                                                        mpasOcean.cellsOnCell,
                                                                        mpasOcean.areaCell,
                                                                        mpasOcean.edgeSignOnCell,
                                                                        mpasOcean.dvEdge,
                                                                        mpasOcean.maxLevelCell)
end

function calculate_thickness_tendency_cuda_kernel!(nCells,
                                              layerThicknessTendency,
                                              layerThickness,
                                              normalVelocity,
                                              maxLevelEdgeTop,
                                              nEdgesOnCell,
                                              edgesOnCell,
                                              cellsOnCell,
                                              areaCell,
                                              edgeSignOnCell,
                                              dvEdge,
                                              maxLevelCell)

    iCell = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
    if iCell <= nCells
        
        for k in 1:maxLevelCell[iCell]
            layerThicknessTendency[k,iCell] = 0
        end
        
        # sum flux through each edge of cell
        for i in 1:nEdgesOnCell[iCell]
            iEdge = edgesOnCell[i,iCell]
            iCell2 = cellsOnCell[i,iCell]
            
            for k in 1:maxLevelEdgeTop[iEdge]
                layerThicknessEdge = 0.5 * ( layerThickness[k,iCell2] + layerThickness[k,iCell] )
                layerThicknessTendency[k,iCell] += normalVelocity[k,iEdge] * edgeSignOnCell[iCell,i] * dvEdge[iEdge] * layerThicknessEdge / areaCell[iCell]
            end
        end
        
    end
    return
end


function update_thickness_by_tendency_cuda!(mpasOcean::MPAS_Ocean_CUDA)
    CUDA.@cuda blocks=cld(mpasOcean.nCells, 1024) threads=1024 maxregs=64 update_thickness_by_tendency_cuda_kernel!(
                                                                        mpasOcean.nCells,
                                                                        mpasOcean.layerThickness,
                                                                        mpasOcean.dt,
                                                                        mpasOcean.layerThicknessTendency,
                                                                        mpasOcean.maxLevelCell)
end

function update_thickness_by_tendency_cuda_kernel!(nCells,
                                             layerThickness,
                                             dt,
                                             layerThicknessTendency,
                                             maxLevelCell)
    iCell = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
    if iCell <= nCells
        
        for k in 1:maxLevelCell[iCell]
            layerThickness[k,iCell] += dt * layerThicknessTendency[k,iCell]
        end

    end
    return
end



