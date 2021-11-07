import CUDA

include("../mode_init/MPAS_Ocean.jl")

### CPU tendency calculation

function calculate_normal_velocity_tendency!(mpasOcean::MPAS_Ocean)
    mpasOcean.normalVelocityTendency[:,:] .= 0

    for iEdge in 1:mpasOcean.nEdges
        # gravity term
        cell1Index, cell2Index = mpasOcean.cellsOnEdge[:,iEdge]
        
        if cell1Index != 0 && cell2Index != 0
            # gradient of ssh across edge
            mpasOcean.normalVelocityTendency[iEdge,:] .= mpasOcean.gravity * ( mpasOcean.sshCurrent[cell1Index] - mpasOcean.sshCurrent[cell2Index] ) / mpasOcean.dcEdge[iEdge]
        end
        
        # coriolis term
        for i in 1:mpasOcean.nEdgesOnEdge[iEdge]
            eoe = mpasOcean.edgesOnEdge[i,iEdge]
            
            if eoe != 0
                for k in 1:mpasOcean.maxLevelEdgeTop[iEdge]+1
                    mpasOcean.normalVelocityTendency[iEdge,k] += mpasOcean.weightsOnEdge[i,iEdge] * mpasOcean.normalVelocityCurrent[eoe,k] * mpasOcean.fEdge[eoe]
                    mpasOcean.normalVelocityTendency[iEdge,k] *= mpasOcean.edgeMask[iEdge,k]
                end
            end
        end
    end
end

function update_normal_velocity_by_tendency!(mpasOcean::MPAS_Ocean)
    mpasOcean.normalVelocityCurrent .+= mpasOcean.dt .* mpasOcean.normalVelocityTendency
end



function calculate_ssh_tendency!(mpasOcean::MPAS_Ocean)
    mpasOcean.sshTendency[:] .= 0
    
    for iCell in 1:mpasOcean.nCells
        # sum flux through each edge of cell
        netflux = 0
        for i in 1:mpasOcean.nEdgesOnCell[iCell]
            edgeID = mpasOcean.edgesOnCell[i,iCell]
            neighborCellID = mpasOcean.cellsOnCell[i,iCell]
            
             # flux in top layer computed seperately because of variable height
            if neighborCellID != 0
               
                sshedge = 0.5*( mpasOcean.sshCurrent[iCell] + mpasOcean.sshCurrent[neighborCellID] )
            else
                sshedge = mpasOcean.sshCurrent[iCell]
            end
            
            flux = ( mpasOcean.layerThicknessEdge[edgeID,1] + sshedge ) * mpasOcean.normalVelocityCurrent[edgeID,1] * mpasOcean.dvEdge[edgeID]
            netflux += flux * mpasOcean.edgeSignOnCell[iCell,i]
            
            # lower layers
            for k in 2:mpasOcean.maxLevelEdgeTop[edgeID]+1
                flux = (mpasOcean.layerThicknessEdge[edgeID,k]) * mpasOcean.normalVelocityCurrent[edgeID,k] * mpasOcean.dvEdge[edgeID]
                netflux += flux * mpasOcean.edgeSignOnCell[iCell,i]
            end
        end
        mpasOcean.sshTendency[iCell] = netflux / mpasOcean.areaCell[iCell]
    end
end

function update_ssh_by_tendency!(mpasOcean::MPAS_Ocean)
    mpasOcean.sshCurrent .+= mpasOcean.dt .* mpasOcean.sshTendency
end





### GPU tendency calculation


function calculate_normal_velocity_tendency_cuda!(mpasOcean::MPAS_Ocean)
    CUDA.@cuda blocks=cld(mpasOcean.nEdges, 1024) threads=1024 maxregs=64 calculate_normal_velocity_tendency_cuda_kernel!(
                                                                        mpasOcean.nEdges,
                                                                        mpasOcean.normalVelocityTendency,
                                                                        mpasOcean.normalVelocityCurrent,
                                                                        mpasOcean.sshCurrent,
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
                                                          ssh,
                                                          cellsOnEdge,
                                                          nEdgesOnEdge,
                                                          edgesOnEdge,
                                                          weightsOnEdge,
                                                          fEdge,
                                                          dcEdge,
                                                          gravity)
    
    iEdge = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
    if iEdge <= nEdges

        
        # gravity term: take gradient of ssh across edge
        normalVelocityTendency[iEdge] = gravity * ( ssh[cellsOnEdge[1,iEdge]] - ssh[cellsOnEdge[2,iEdge]] ) / dcEdge[iEdge]

        # coriolis term: combine norm. vel. of surrounding edges to approx. tangential vel.
        for i = 1:nEdgesOnEdge[iEdge]
            normalVelocityTendency[iEdge] += weightsOnEdge[i,iEdge] * normalVelocity[edgesOnEdge[i,iEdge]] * fEdge[edgesOnEdge[i,iEdge]]
        end

        
    end
    return
end







function update_normal_velocity_by_tendency_cuda!(mpasOcean::MPAS_Ocean)
    CUDA.@cuda blocks=cld(mpasOcean.nEdges, 1024) threads=1024 maxregs=64 update_normal_velocity_by_tendency_cuda_kernel!(
                                                                        mpasOcean.nEdges,
                                                                        mpasOcean.normalVelocityCurrent,
                                                                        mpasOcean.dt,
                                                                        mpasOcean.normalVelocityTendency)
end

function update_normal_velocity_by_tendency_cuda_kernel!(nEdges,
                                                         normalVelocityCurrent,
                                                         dt,
                                                         normalVelocityTendency)
    iEdge = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
    if iEdge <= nEdges
        
        normalVelocityCurrent[iEdge] += dt * normalVelocityTendency[iEdge]
        
    end
    return
end
    







function calculate_ssh_tendency_cuda!(mpasOcean::MPAS_Ocean)
    CUDA.@cuda blocks=cld(mpasOcean.nCells, 1024) threads=1024 maxregs=64 calculate_ssh_tendency_cuda_kernel!(
                                                                        mpasOcean.nCells,
                                                                        mpasOcean.sshTendency,
                                                                        mpasOcean.sshCurrent,
                                                                        mpasOcean.normalVelocityCurrent,
                                                                        mpasOcean.bottomDepth,
                                                                        mpasOcean.nEdgesOnCell,
                                                                        mpasOcean.edgesOnCell,
                                                                        mpasOcean.cellsOnCell,
                                                                        mpasOcean.areaCell,
                                                                        mpasOcean.edgeSignOnCell,
                                                                        mpasOcean.dvEdge)
end

function calculate_ssh_tendency_cuda_kernel!(nCells,
                                              sshTendency,
                                              ssh,
                                              normalVelocity,
                                              bottomDepth,
                                              nEdgesOnCell,
                                              edgesOnCell,
                                              cellsOnCell,
                                              areaCell,
                                              edgeSignOnCell,
                                              dvEdge)
    
    iCell = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
    if iCell <= nCells
        
        sshTendency[iCell] = 0
        
        # sum flux through each edge of cell
        for i in 1:nEdgesOnCell[iCell]
            edgeID = edgesOnCell[i,iCell]
            neighborCellID = cellsOnCell[i,iCell]

            if neighborCellID !== 0
                depth = ( bottomDepth[neighborCellID] + bottomDepth[iCell] ) / 2.0
            else
                depth = bottomDepth[iCell]
            end
            
            flux = depth * normalVelocity[edgeID]
            sshTendency[iCell] += flux * edgeSignOnCell[iCell,i] * dvEdge[edgeID] / areaCell[iCell]
        end
        
    end
    return
end







function update_ssh_by_tendency_cuda!(mpasOcean::MPAS_Ocean)
    CUDA.@cuda blocks=cld(mpasOcean.nCells, 1024) threads=1024 maxregs=64 update_ssh_by_tendency_cuda_kernel!(
                                                                        mpasOcean.nCells,
                                                                        mpasOcean.sshCurrent,
                                                                        mpasOcean.dt,
                                                                        mpasOcean.sshTendency)
end

function update_ssh_by_tendency_cuda_kernel!(nCells,
                                             sshCurrent,
                                             dt,
                                             sshTendency)
    iCell = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
    if iCell <= nCells
        
        sshCurrent[iCell] += dt * sshTendency[iCell]
        
    end
    return
end








