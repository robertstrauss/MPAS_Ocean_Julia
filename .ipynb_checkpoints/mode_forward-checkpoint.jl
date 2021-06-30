using CUDA
using KernelAbstractions


function ComputeSSHTendency_singleindex(iCell::Integer, myMPAS_O, normalVelocity, ssh)
    ## dn/dt = - divergence(h*u)
    SSHTendency = 0

    # sum flux through each edge of cell
    for i in 1:myMPAS_O.nEdgesOnCell[iCell]
        edgeID = myMPAS_O.edgesOnCell[i,iCell]
        neighborCellID = myMPAS_O.cellsOnCell[i,iCell]

        if neighborCellID !== 0 && edgeID !== 0 # boundary conditions: zero flux at boundaries
            mean_depth = ( myMPAS_O.bottomDepth[neighborCellID] + myMPAS_O.bottomDepth[iCell] ) / 2.0
            flux = mean_depth * normalVelocity[edgeID]
            SSHTendency += flux * myMPAS_O.edgeSignOnCell[iCell,i] * myMPAS_O.dvEdge[edgeID]
        end
    end
    SSHTendency /= myMPAS_O.areaCell[iCell]

    return SSHTendency
end

function ComputeSSHTendency_func(myMPAS_O::MPAS_Ocean, normalVelocity, ssh)
    ## dn/dt = - divergence(h*u)
    SSHTendency = zeros(Float64, myMPAS_O.nCells)

    for iCell in 1:myMPAS_O.nCells
        SSHTendency[iCell] = ComputeSSHTendency(iCell, myMPAS_O, normalVelocity, ssh)
    end

    return SSHTendency
end


function ComputeSSHTendency(myMPAS_O::MPAS_Ocean, normalVelocity, ssh)
    ## dn/dt = - divergence(h*u)
    SSHTendency = zeros(Float64, myMPAS_O.nCells)

    for iCell in 1:myMPAS_O.nCells
        # sum flux through each edge of cell
        for i in 1:myMPAS_O.nEdgesOnCell[iCell]
            edgeID = myMPAS_O.edgesOnCell[i,iCell]
            neighborCellID = myMPAS_O.cellsOnCell[i,iCell]
            
            if neighborCellID !== 0
                mean_depth = ( myMPAS_O.bottomDepth[neighborCellID] + myMPAS_O.bottomDepth[iCell] ) / 2.0
            else
                mean_depth = myMPAS_O.bottomDepth[iCell]
            end
            
            flux = mean_depth * normalVelocity[edgeID]
            SSHTendency[iCell] += flux * myMPAS_O.edgeSignOnCell[iCell,i] * myMPAS_O.dvEdge[edgeID]
        end
        SSHTendency[iCell] /= myMPAS_O.areaCell[iCell]
    end

    return SSHTendency
end







function ComputeNormalVelocityTendency_singleindex(iEdge::Integer, myMPAS_O, normalVelocity, ssh)
    ## du/dt = - g dn/dx + coriolis + other
    normalVelocityTendency = 0
    
    # take gradient of sshCurrent across edge
    cell1Index, cell2Index = myMPAS_O.cellsOnEdge[:,iEdge]

    if cell1Index !== 0 && cell2Index !== 0
        # TODO add nonlinear terms
        normalVelocityTendency[iEdge] = gravity * ( ssh[cell1Index] - ssh[cell2Index] ) / myMPAS_O.dcEdge[iEdge] # gravity term
    else
        # TODO boundary condition
    end

    # coriolis term
    for i in 1:myMPAS_O.nEdgesOnEdge[iEdge]
        eoe = myMPAS_O.edgesOnEdge[i,iEdge]

        if eoe !== 0 # boundary conditions: coriolis not influenced by edges outside boundaries
            normalVelocityTendency += myMPAS_O.weightsOnEdge[i,iEdge] * normalVelocity[eoe] * myMPAS_O.fEdge[eoe]
        end
    end

    return normalVelocityTendency
end


function ComputeNormalVelocityTendency_func(myMPAS_O::MPAS_Ocean, normalVelocity, ssh)
    normalVelocityTendency = zeros(myMPAS_O.nEdges)
    for iEdge in 1:myMPAS_O.nEdges
        normalVelocityTendency = ComputeNormalVelocityTendency(iEdge, myMPAS_O, normalVelocity, ssh)
    end
    return normalVelocityTendency
end

function ComputeNormalVelocityTendency(myMPAS_O::MPAS_Ocean, normalVelocity, ssh)
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













function forwardEuler!(myMPAS_O::MPAS_Ocean, dt, time)
    for t in 0:dt:time
        myMPAS_O.normalVelocityNew = myMPAS_O.normalVelocityCurrent + dt*ComputeNormalVelocityTendency(myMPAS_O, myMPAS_O.normalVelocityCurrent, myMPAS_O.sshCurrent)
        myMPAS_O.sshNew = myMPAS_O.sshCurrent + dt*ComputeSSHTendency(myMPAS_O, myMPAS_O.normalVelocityCurrent, myMPAS_O.sshCurrent)

        myMPAS_O.normalVelocityCurrent[:] = myMPAS_O.normalVelocityNew[:]
        myMPAS_O.sshCurrent[:] = myMPAS_O.sshNew[:]
    end
end


function forwardBackward!(myMPAS_O::MPAS_Ocean, dt, time)
    for t in 0:dt:time
        myMPAS_O.normalVelocityNew = myMPAS_O.normalVelocityCurrent + dt*ComputeNormalVelocityTendency(myMPAS_O, myMPAS_O.normalVelocityCurrent, myMPAS_O.sshCurrent)
        myMPAS_O.sshNew = myMPAS_O.sshCurrent + dt*ComputeSSHTendency(myMPAS_O, myMPAS_O.normalVelocityNew, myMPAS_O.sshCurrent)
        
        
        
        myMPAS_O.normalVelocityCurrent[:] = myMPAS_O.normalVelocityNew[:]
        myMPAS_O.sshCurrent[:] = myMPAS_O.sshNew[:]
        

    end
end


function forwardBackwardSingleStep!(myMPAS_O, dt)
    myMPAS_O.normalVelocityNew = myMPAS_O.normalVelocityCurrent + dt*ComputeNormalVelocityTendency(myMPAS_O, myMPAS_O.normalVelocityCurrent, myMPAS_O.sshCurrent)
    myMPAS_O.sshNew = myMPAS_O.sshCurrent + dt*ComputeSSHTendency(myMPAS_O, myMPAS_O.normalVelocityNew, myMPAS_O.sshCurrent)

    myMPAS_O.normalVelocityCurrent[:] = myMPAS_O.normalVelocityNew[:]
    myMPAS_O.sshCurrent[:] = myMPAS_O.sshNew[:]
end




@kernel function stepNormalVelocity!(normalVelocityNew, normalVelocity, ssh, cellsOnEdge, nEdgesOnEdge, edgesOnEdge, weightsOnEdge, fEdge, dcEdge, dt)
    
    
    
    iEdge = @index(Global)
    ## du/dt = - g dn/dx + coriolis + other
    
    
    
    normalVelocityNew[iEdge] = normalVelocity[iEdge]
    
    # take gradient of sshCurrent across edge
    cell1Index, cell2Index = cellsOnEdge[:,iEdge]

    if cell1Index !== 0 && cell2Index !== 0
        # TODO add nonlinear terms
        normalVelocityNew[iEdge] += gravity * ( ssh[cell1Index] - ssh[cell2Index] ) / dcEdge[iEdge] # gravity term
    else
        # TODO boundary condition
    end
    
    
    

    # coriolis term
    for i in 1:nEdgesOnEdge[iEdge]
        eoe = edgesOnEdge[i,iEdge]

        if eoe !== 0 # boundary conditions: coriolis not influenced by edges outside boundaries
            normalVelocityNew[iEdge] += weightsOnEdge[i,iEdge] * normalVelocity[eoe] * fEdge[eoe]
        end
    end
    
#     @inbounds myMPAS_O.normalVelocityNew[iEdge] = myMPAS_O.normalVelocityCurrent[iEdge] + dt* normalVelocityTendency#ComputeNormalVelocityTendency(iEdge, myMPAS_O, normalVelocity, ssh)
    
    nothing
end
    
@kernel function stepSSH!(dt)
#     iCell = @index(Global)
#     @inbounds myMPAS_O.sshNew[iEdge] = myMPAS_O.sshCurrent[iEdge] + dt*ComputeSSHTendency(iCell, myMPAS_O, normalVelocity, ssh)
    
    nothing
end


function forwardBackwardStep!(myMPAS_O, dt)
    
#     normalVelocity = CuArray(myMPAS_O.normalVelocityCurrent)
#     ssh = CuArray(myMPAS_O.sshCurrent)
#     nEdgesOnCell = CuArray(myMPAS_O.nEdgesOnCell)
#     edgesOnCell = CuArray(myMPAS_O.edgesOnCell)
#     cellsOnCell = CuArray(myMPAS_O.cellsOnCell)
#     bottomDepth = CuArray(myMPAS_O.bottomDepth)
#     edgeSignOnCell = CuArray(myMPAS_O.edgeSignOnCell)
#     dvEdge = CuArray(myMPAS_O.dvEdge)
#     areaCell = CuArray(myMPAS_O.areaCell)
    
    stepNormalVelocityKernel! = stepNormalVelocity!(myMPAS_O.device, myMPAS_O.workGroupSize)
    stepSSHKernel!            = stepSSH!(myMPAS_O.device, myMPAS_O.workGroupSize)
    
    stepNormalVelocityEvent = stepNormalVelocityKernel!(myMPAS_O.normalVelocityNew,
                                                        myMPAS_O.normalVelocityCurrent,
                                                        myMPAS_O.sshCurrent,
                                                        myMPAS_O.cellsOnEdge,
                                                        myMPAS_O.nEdgesOnEdge,
                                                        myMPAS_O.edgesOnEdge,
                                                        myMPAS_O.weightsOnEdge,
                                                        myMPAS_O.fEdge,
                                                        myMPAS_O.dcEdge,
                                                        dt, ndrange=myMPAS_O.nEdges)
    wait(myMPAS_O.device, stepNormalVelocityEvent)
    
    stepSSHEvent            = stepSSHKernel!(dt, ndrange=myMPAS_O.nCells)
    wait(myMPAS_O.device, stepSSHEvent)
    
    return nothing
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
        
        