


function ComputeSSHTendency(myMPAS_O, normalVelocity, ssh)
    ## dn/dt = - divergence(h*u)
    SSHTendency = zeros(Float64, myMPAS_O.nCells)

    for iCell in 1:myMPAS_O.nCells
        # sum flux through each edge of cell
        for i in 1:myMPAS_O.nEdgesOnCell[iCell]
            edgeID = myMPAS_O.edgesOnCell[i,iCell]
            neighborCellID = myMPAS_O.cellsOnCell[i,iCell]
            
            if neighborCellID !== 0 && edgeID !== 0 # boundary conditions: zero flux at boundaries
                mean_depth = ( myMPAS_O.bottomDepth[neighborCellID] + myMPAS_O.bottomDepth[iCell] ) / 2.0
                flux = mean_depth * normalVelocity[edgeID] # TODO multiply by edge length
                SSHTendency[iCell] += flux * myMPAS_O.edgeSignOnCell[iCell,i] * myMPAS_O.dvEdge[edgeID]
            end
        end
        SSHTendency[iCell] /= myMPAS_O.areaCell[iCell]
    end

    return SSHTendency
end



function ComputeNormalVelocityTendency(myMPAS_O, normalVelocity, ssh)
    ## du/dt = - g dn/dx + coriolis + other
    normalVelocityTendency = zeros(Float64, myMPAS_O.nEdges)

    for iEdge in 1:myMPAS_O.nEdges
        # take gradient of sshCurrent across edge
        cell1Index, cell2Index = myMPAS_O.cellsOnEdge[:,iEdge]
        
        if cell1Index !== 0 && cell2Index !== 0
            # TODO add nonlinear terms
            normalVelocityTendency[iEdge] = gravity * ( ssh[cell1Index] - ssh[cell2Index] ) / myMPAS_O.dcEdge[iEdge] # gravity term
        end
        
        # coriolis term
        for i in 1:myMPAS_O.nEdgesOnEdge[iEdge]
            eoe = myMPAS_O.edgesOnEdge[i,iEdge]
            
            if eoe !== 0 # boundary conditions: coriolis not influenced by edges outside boundaries
                normalVelocityTendency[iEdge] += myMPAS_O.weightsOnEdge[i,iEdge] * normalVelocity[eoe] * myMPAS_O.fEdge[eoe]
            end
        end
    end

    return normalVelocityTendency
end



function forwardEuler!(myMPAS_O, dt, time)
    for t in 0:dt:time
        myMPAS_O.normalVelocityNew = myMPAS_O.normalVelocityCurrent + dt*ComputeNormalVelocityTendency(myMPAS_O, myMPAS_O.normalVelocityCurrent, myMPAS_O.sshCurrent)
        myMPAS_O.sshNew = myMPAS_O.sshCurrent + dt*ComputeSSHTendency(myMPAS_O, myMPAS_O.normalVelocityCurrent, myMPAS_O.sshCurrent)

        myMPAS_O.normalVelocityCurrent[:] = myMPAS_O.normalVelocityNew[:]
        myMPAS_O.sshCurrent[:] = myMPAS_O.sshNew[:]
    end
end


function forwardBackward!(myMPAS_O, dt, time)
    for t in 0:dt:time
        myMPAS_O.normalVelocityNew = myMPAS_O.normalVelocityCurrent + dt*ComputeNormalVelocityTendency(myMPAS_O, myMPAS_O.normalVelocityCurrent, myMPAS_O.sshCurrent)
        myMPAS_O.sshNew = myMPAS_O.sshCurrent + dt*ComputeSSHTendency(myMPAS_O, myMPAS_O.normalVelocityNew, myMPAS_O.sshCurrent)
        
        
        
        myMPAS_O.normalVelocityCurrent[:] = myMPAS_O.normalVelocityNew[:]
        myMPAS_O.sshCurrent[:] = myMPAS_O.sshNew[:]
        
        boundaryConditions!(myMPAS_O)
    end
end


# function integrate!(method, myMPAS_O, dt, simulationTime; callback = function() end, callbackTimes=[])
#     discreteCallbackTimes = round.(collect(callbackTimes)/dt)*dt
    
#     simulationDurations = discreteCallbackTimes[2:end] - discreteCallbackTimes[1:end-1]
#     append!(simulationDurations, 0)
    
#     for i in 1:length(discreteCallbackTimes)
#         callback(myMPAS_O, i, discreteCallbackTimes[i])
#         method(myMPAS_O, dt, simulationDurations[i])
#     end
    
# end
        
        