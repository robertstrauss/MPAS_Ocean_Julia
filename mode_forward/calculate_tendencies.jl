include("../mode_init/MPAS_Ocean.jl")

### CPU tendency calculation

function calculate_normal_velocity_tendency!(mpasOcean::MPAS_Ocean)
    @fastmath for iEdge in 1:mpasOcean.nEdges
		mpasOcean.normalVelocityTendency[iEdge,:] .= 0.0

        if mpasOcean.boundaryEdge[iEdge] == 0
            # gravity term: take gradient of sshCurrent across edge
            cell1Index = mpasOcean.cellsOnEdge[1,iEdge]
			cell2Index = mpasOcean.cellsOnEdge[2,iEdge]

            if cell1Index != 0 && cell2Index != 0
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
	@fastmath for iEdge::Int64 in 1:mpasOcean.nEdges
		mpasOcean.normalVelocityCurrent[iEdge] += mpasOcean.dt * mpasOcean.normalVelocityTendency[iEdge]
	end
end

function calculate_ssh_tendency!(mpasOcean::MPAS_Ocean)
    @fastmath for iCell in 1:mpasOcean.nCells # Threads.@threads
    	mpasOcean.sshTendency[iCell] = 0.0
        # sum flux through each edge of cell
        netflux = 0.0
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
