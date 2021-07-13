using CUDA
# using KernelAbstractions


function ComputeSSHTendency_singleindex(iCell::Integer, mpasOcean, normalVelocity, ssh)
    ## dn/dt = - divergence(h*u)
    SSHTendency = 0

    # sum flux through each edge of cell
    for i in 1:mpasOcean.nEdgesOnCell[iCell]
        edgeID = mpasOcean.edgesOnCell[i,iCell]
        neighborCellID = mpasOcean.cellsOnCell[i,iCell]

        if neighborCellID !== 0 && edgeID !== 0 # boundary conditions: zero flux at boundaries
            mean_depth = ( mpasOcean.bottomDepth[neighborCellID] + mpasOcean.bottomDepth[iCell] ) / 2.0
            flux = mean_depth * normalVelocity[edgeID]
            SSHTendency += flux * mpasOcean.edgeSignOnCell[iCell,i] * mpasOcean.dvEdge[edgeID]
        end
    end
    SSHTendency /= mpasOcean.areaCell[iCell]

    return SSHTendency
end

function ComputeSSHTendency_func(mpasOcean::MPAS_Ocean, normalVelocity, ssh)
    ## dn/dt = - divergence(h*u)
    SSHTendency = zeros(Float64, mpasOcean.nCells)

    for iCell in 1:mpasOcean.nCells
        SSHTendency[iCell] = ComputeSSHTendency(iCell, mpasOcean, normalVelocity, ssh)
    end

    return SSHTendency
end


function ComputeSSHTendency(mpasOcean::MPAS_Ocean, normalVelocity, ssh)
    ## dn/dt = - divergence(h*u)
    SSHTendency = zeros(Float64, mpasOcean.nCells)

    for iCell in 1:mpasOcean.nCells
        # sum flux through each edge of cell
        for i in 1:mpasOcean.nEdgesOnCell[iCell]
            edgeID = mpasOcean.edgesOnCell[i,iCell]
            neighborCellID = mpasOcean.cellsOnCell[i,iCell]
            
            if neighborCellID !== 0
                mean_depth = ( mpasOcean.bottomDepth[neighborCellID] + mpasOcean.bottomDepth[iCell] ) / 2.0
            else
                mean_depth = mpasOcean.bottomDepth[iCell]
            end
            
            flux = mean_depth * normalVelocity[edgeID]
            SSHTendency[iCell] += flux * mpasOcean.edgeSignOnCell[iCell,i] * mpasOcean.dvEdge[edgeID]
        end
        SSHTendency[iCell] /= mpasOcean.areaCell[iCell]
    end

    return SSHTendency
end







function ComputeNormalVelocityTendency_singleindex(iEdge::Integer, mpasOcean, normalVelocity, ssh)
    ## du/dt = - g dn/dx + coriolis + other
    normalVelocityTendency = 0
    
    # take gradient of sshCurrent across edge
    cell1Index, cell2Index = mpasOcean.cellsOnEdge[:,iEdge]

    if cell1Index !== 0 && cell2Index !== 0
        normalVelocityTendency[iEdge] = gravity * ( ssh[cell1Index] - ssh[cell2Index] ) / mpasOcean.dcEdge[iEdge] # gravity term
    end

    # coriolis term
    for i in 1:mpasOcean.nEdgesOnEdge[iEdge]
        eoe = mpasOcean.edgesOnEdge[i,iEdge]

        if eoe !== 0 # boundary conditions: coriolis not influenced by edges outside boundaries
            normalVelocityTendency += mpasOcean.weightsOnEdge[i,iEdge] * normalVelocity[eoe] * mpasOcean.fEdge[eoe]
        end
    end

    return normalVelocityTendency
end


function ComputeNormalVelocityTendency_func(mpasOcean::MPAS_Ocean, normalVelocity, ssh)
    normalVelocityTendency = zeros(mpasOcean.nEdges)
    for iEdge in 1:mpasOcean.nEdges
        normalVelocityTendency = ComputeNormalVelocityTendency(iEdge, mpasOcean, normalVelocity, ssh)
    end
    return normalVelocityTendency
end

function ComputeNormalVelocityTendency(mpasOcean::MPAS_Ocean, normalVelocity, ssh)
    ## du/dt = - g dn/dx + coriolis + other
    normalVelocityTendency = zeros(Float64, mpasOcean.nEdges)

    for iEdge in 1:mpasOcean.nEdges
        # take gradient of sshCurrent across edge
        cell1Index, cell2Index = mpasOcean.cellsOnEdge[:,iEdge]
        
        if cell1Index !== 0 && cell2Index !== 0
        normalVelocityTendency[iEdge] = gravity * ( ssh[cell1Index] - ssh[cell2Index] ) / mpasOcean.dcEdge[iEdge] # gravity term
        end
        
        # coriolis term
        for i in 1:mpasOcean.nEdgesOnEdge[iEdge]
            eoe = mpasOcean.edgesOnEdge[i,iEdge]
            
            if eoe !== 0 # boundary conditions: coriolis not influenced by edges outside boundaries
                normalVelocityTendency[iEdge] += mpasOcean.weightsOnEdge[i,iEdge] * normalVelocity[eoe] * mpasOcean.fEdge[eoe]
            end
        end
    end

    return normalVelocityTendency
end













function forwardEuler!(mpasOcean::MPAS_Ocean, dt, time)
    for t in 0:dt:time
        mpasOcean.normalVelocityNew = mpasOcean.normalVelocityCurrent + dt*ComputeNormalVelocityTendency(mpasOcean, mpasOcean.normalVelocityCurrent, mpasOcean.sshCurrent)
        mpasOcean.sshNew = mpasOcean.sshCurrent + dt*ComputeSSHTendency(mpasOcean, mpasOcean.normalVelocityCurrent, mpasOcean.sshCurrent)

        mpasOcean.normalVelocityCurrent[:] = mpasOcean.normalVelocityNew[:]
        mpasOcean.sshCurrent[:] = mpasOcean.sshNew[:]
    end
end


function forwardBackward!(mpasOcean::MPAS_Ocean, dt, time)
    for t in 0:dt:time
        mpasOcean.normalVelocityNew = mpasOcean.normalVelocityCurrent + dt*ComputeNormalVelocityTendency(mpasOcean, mpasOcean.normalVelocityCurrent, mpasOcean.sshCurrent)
        mpasOcean.sshNew = mpasOcean.sshCurrent + dt*ComputeSSHTendency(mpasOcean, mpasOcean.normalVelocityNew, mpasOcean.sshCurrent)
        
        
        
        mpasOcean.normalVelocityCurrent[:] = mpasOcean.normalVelocityNew[:]
        mpasOcean.sshCurrent[:] = mpasOcean.sshNew[:]
        

    end
end


function forwardBackwardSingleStep!(mpasOcean, dt)
    mpasOcean.normalVelocityNew = mpasOcean.normalVelocityCurrent + dt*ComputeNormalVelocityTendency(mpasOcean, mpasOcean.normalVelocityCurrent, mpasOcean.sshCurrent)
    mpasOcean.sshNew = mpasOcean.sshCurrent + dt*ComputeSSHTendency(mpasOcean, mpasOcean.normalVelocityNew, mpasOcean.sshCurrent)

    mpasOcean.normalVelocityCurrent[:] = mpasOcean.normalVelocityNew[:]
    mpasOcean.sshCurrent[:] = mpasOcean.sshNew[:]
end




# @kernel function stepNormalVelocity!(normalVelocityNew, normalVelocity, ssh, cellsOnEdge, nEdgesOnEdge, edgesOnEdge, weightsOnEdge, fEdge, dcEdge, dt)
    
    
    
#     iEdge = @index(Global)
#     ## du/dt = - g dn/dx + coriolis + other
    
    
    
#     normalVelocityNew[iEdge] = normalVelocity[iEdge]
    
#     # take gradient of sshCurrent across edge
#     cell1Index, cell2Index = cellsOnEdge[:,iEdge]

#     if cell1Index !== 0 && cell2Index !== 0
#         normalVelocityNew[iEdge] += gravity * ( ssh[cell1Index] - ssh[cell2Index] ) / dcEdge[iEdge] # gravity term
#     end
    
    
    

#     # coriolis term
#     for i in 1:nEdgesOnEdge[iEdge]
#         eoe = edgesOnEdge[i,iEdge]

#         if eoe !== 0 # boundary conditions: coriolis not influenced by edges outside boundaries
#             normalVelocityNew[iEdge] += weightsOnEdge[i,iEdge] * normalVelocity[eoe] * fEdge[eoe]
#         end
#     end
    
# #     @inbounds mpasOcean.normalVelocityNew[iEdge] = mpasOcean.normalVelocityCurrent[iEdge] + dt* normalVelocityTendency#ComputeNormalVelocityTendency(iEdge, mpasOcean, normalVelocity, ssh)
    
#     nothing
# end
    
# @kernel function stepSSH!(dt)
# #     iCell = @index(Global)
# #     @inbounds mpasOcean.sshNew[iEdge] = mpasOcean.sshCurrent[iEdge] + dt*ComputeSSHTendency(iCell, mpasOcean, normalVelocity, ssh)
    
#     nothing
# end


# function forwardBackwardStep!(mpasOcean, dt)
    
# #     normalVelocity = CuArray(mpasOcean.normalVelocityCurrent)
# #     ssh = CuArray(mpasOcean.sshCurrent)
# #     nEdgesOnCell = CuArray(mpasOcean.nEdgesOnCell)
# #     edgesOnCell = CuArray(mpasOcean.edgesOnCell)
# #     cellsOnCell = CuArray(mpasOcean.cellsOnCell)
# #     bottomDepth = CuArray(mpasOcean.bottomDepth)
# #     edgeSignOnCell = CuArray(mpasOcean.edgeSignOnCell)
# #     dvEdge = CuArray(mpasOcean.dvEdge)
# #     areaCell = CuArray(mpasOcean.areaCell)
    
#     stepNormalVelocityKernel! = stepNormalVelocity!(mpasOcean.device, mpasOcean.workGroupSize)
#     stepSSHKernel!            = stepSSH!(mpasOcean.device, mpasOcean.workGroupSize)
    
#     stepNormalVelocityEvent = stepNormalVelocityKernel!(mpasOcean.normalVelocityNew,
#                                                         mpasOcean.normalVelocityCurrent,
#                                                         mpasOcean.sshCurrent,
#                                                         mpasOcean.cellsOnEdge,
#                                                         mpasOcean.nEdgesOnEdge,
#                                                         mpasOcean.edgesOnEdge,
#                                                         mpasOcean.weightsOnEdge,
#                                                         mpasOcean.fEdge,
#                                                         mpasOcean.dcEdge,
#                                                         dt, ndrange=mpasOcean.nEdges)
#     wait(mpasOcean.device, stepNormalVelocityEvent)
    
#     stepSSHEvent            = stepSSHKernel!(dt, ndrange=mpasOcean.nCells)
#     wait(mpasOcean.device, stepSSHEvent)
    
#     return nothing
# end

# function integrate!(method, mpasOcean, dt, simulationTime; callback = function() end, callbackTimes=[])
#     discreteCallbackTimes = round.(collect(callbackTimes)/dt)*dt
    
#     simulationDurations = discreteCallbackTimes[2:end] - discreteCallbackTimes[1:end-1]
#     append!(simulationDurations, 0)
    
#     for i in 1:length(discreteCallbackTimes)
#         callback(mpasOcean, i, discreteCallbackTimes[i])
#         method(mpasOcean, dt, simulationDurations[i])
#     end
    
# end
        
        


# function makeStepKernels(mpasOcean)
#     stepNormalVelocityKernel! = stepNormalVelocity!(mpasOcean.device, mpasOcean.workGroupSize)
#     stepSSHKernel!            = stepSSH!(mpasOcean.device, mpasOcean.workGroupSize)
#     return stepNormalVelocityKernel!, stepSSHKernel!
# end

# function forwardBackwardStep!(mpasOcean, dt, stepNormalVelocityKernel!, stepSSHKernel!)
#     stepNormalVelocityEvent = stepNormalVelocityKernel!(mpasOcean.normalVelocityNew,
#                                                         mpasOcean.normalVelocityCurrent,
#                                                         mpasOcean.sshCurrent,
#                                                         mpasOcean.cellsOnEdge,
#                                                         mpasOcean.nEdgesOnEdge,
#                                                         mpasOcean.edgesOnEdge,
#                                                         mpasOcean.weightsOnEdge,
#                                                         mpasOcean.fEdge,
#                                                         mpasOcean.dcEdge,
#                                                         dt, ndrange=mpasOcean.nEdges)
#     wait(mpasOcean.device, stepNormalVelocityEvent)

#     stepSSHEvent            = stepSSHKernel!(mpasOcean.sshNew,
#                                              mpasOcean.sshCurrent,
#                                              mpasOcean.normalVelocityNew,
#                                              mpasOcean.bottomDepth,
#                                              mpasOcean.nEdgesOnCell,
#                                              mpasOcean.edgesOnCell,
#                                              mpasOcean.cellsOnCell,
#                                              mpasOcean.areaCell,
#                                              mpasOcean.edgeSignOnCell,
#                                              mpasOcean.dvEdge,
#                                              dt, ndrange=mpasOcean.nCells)
#     wait(mpasOcean.device, stepSSHEvent)
# end