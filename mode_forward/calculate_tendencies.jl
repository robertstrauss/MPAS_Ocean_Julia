import CUDA
import Base.Threads

include("../mode_init/MPAS_Ocean.jl")

### CPU tendency calculation

function calculate_normal_velocity_tendency!(mpasOcean::MPAS_Ocean)
    @inbounds @fastmath Threads.@threads for iEdge::Int64 in 1:mpasOcean.nEdges
		mpasOcean.normalVelocityTendency[iEdge,:] .= 0

        if mpasOcean.boundaryEdge[iEdge] == 0
            # gravity term: take gradient of sshCurrent across edge
            cell1Index = mpasOcean.cellsOnEdge[1,iEdge]
			cell2Index = mpasOcean.cellsOnEdge[2,iEdge]

            if cell1Index != 0 && cell2Index != 0
                mpasOcean.normalVelocityTendency[iEdge,:] .= mpasOcean.gravity * ( mpasOcean.sshCurrent[cell1Index] - mpasOcean.sshCurrent[cell2Index] ) / mpasOcean.dcEdge[iEdge]
            end

            # coriolis term
            for i::Int64 in 1:mpasOcean.nEdgesOnEdge[iEdge]
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
end

function update_normal_velocity_by_tendency_threads!(mpasOcean::MPAS_Ocean)
    Threads.@threads for iEdge in 1:mpasOcean.nEdges
	    mpasOcean.normalVelocityCurrent[iEdge] += mpasOcean.dt * mpasOcean.normalVelocityTendency[iEdge]
    end
end

function update_normal_velocity_by_tendency!(mpasOcean::MPAS_Ocean)
    mpasOcean.normalVelocityCurrent .+= mpasOcean.dt .* mpasOcean.normalVelocityTendency
end

function update_normal_velocity_by_tendency_loop!(mpasOcean::MPAS_Ocean)
	@inbounds @fastmath for iEdge::Int64 in 1:mpasOcean.nEdges
		mpasOcean.normalVelocityCurrent[iEdge] += mpasOcean.dt * mpasOcean.normalVelocityTendency[iEdge]
	end
end

function calculate_ssh_tendency!(mpasOcean::MPAS_Ocean)
    @inbounds @fastmath Threads.@threads for iCell in 1:mpasOcean.nCells
    	mpasOcean.sshTendency[iCell] = 0
        # sum flux through each edge of cell
        netflux = 0
        for i::Int64 in 1:mpasOcean.nEdgesOnCell[iCell]
            edgeID = mpasOcean.edgesOnCell[i,iCell]
            neighborCellID = mpasOcean.cellsOnCell[i,iCell]

            if neighborCellID != 0
                sshedge = 0.5 * ( mpasOcean.sshCurrent[neighborCellID] + mpasOcean.sshCurrent[iCell] )
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



function update_ssh_by_tendency_threads!(mpasOcean::MPAS_Ocean)
    Threads.@threads for iCell in 1:mpasOcean.nCells
    	mpasOcean.sshCurrent .+= mpasOcean.dt .* mpasOcean.sshTendency
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
        cell1 = cellsOnEdge[1,iEdge]
        cell2 = cellsOnEdge[2,iEdge]

        if cell1 != 0 && cell2 != 0
            normalVelocityTendency[iEdge] = gravity * ( ssh[cellsOnEdge[1,iEdge]] - ssh[cellsOnEdge[2,iEdge]] ) / dcEdge[iEdge]
        end

        # coriolis term: combine norm. vel. of surrounding edges to approx. tangential vel.
        for i = 1:nEdgesOnEdge[iEdge]
            eoe = edgesOnEdge[i, iEdge]

            if eoe != 0
                normalVelocityTendency[iEdge] += weightsOnEdge[i,iEdge] * normalVelocity[eoe] * fEdge[eoe]
            end
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

            if neighborCellID != 0
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
