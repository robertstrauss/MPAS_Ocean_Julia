using CUDA

# using KernelAbstractions

using CUDAKernels

CUDA.versioninfo()

CUDA.memory_status()

GPU() = CUDAKernels.CUDADevice()









CODE_ROOT = pwd() * "/"

include(CODE_ROOT * "mode_init.jl")

include(CODE_ROOT * "visualization.jl")

using PyPlot, PyCall
animation  = pyimport("matplotlib.animation")
ipydisplay = pyimport("IPython.display")




# fields = []
# values = []
#
# for field in fieldnames(MPAS_Ocean)
#     if ( typeof(getfield(myMPAS_O, field)) <: CuArray )
#         append!(fields, [field])
#         append!(values, [getfield(myMPAS_O, field)])
#     end
# #     @show typeof(getfield(myMPAS_O, field))
# end
#
#
# append!(fields, [:device])
# append!(values, [myMPAS_O.device])
#
# append!(fields, [:workGroupSize])
# append!(values, [myMPAS_O.workGroupSize])
#
# append!(fields, [:nEdges])
# append!(values, [myMPAS_O.nEdges])
#
# append!(fields, [:nCells])
# append!(values, [myMPAS_O.nCells])
#
# fields = Tuple(fields)
# length(fields), length(values)



function stepNormalVelocityParallel!(normalVelocityNew, normalVelocity, ssh, cellsOnEdge, nEdgesOnEdge, edgesOnEdge, weightsOnEdge, fEdge, dcEdge, dt)
    normalVelocityNew[:] .= 0 #normalVelocity[:]

    # take gradient of sshCurrent across edge
    normalVelocityNew[:] .+= gravity .* ( ssh[cellsOnEdge[1,:]] .- ssh[cellsOnEdge[2,:]] ) ./ dcEdge
    # TODO boundary condition, will error when cell1Index == 0 or cell2Index == 0

    normalVelocity2 = map(iEdge -> iEdge != 0 ? normalVelocity[iEdge] : 0, edgesOnEdge)
    fEdge2 = map(iEdge -> iEdge != 0 ? fEdge[iEdge] : 0, edgesOnEdge)

    # coriolis term

    # edgesOnEdge is a 2x2 matrix, each value is a valid edge index
    # so any edge-shaped thing indexed with it should turn to its shape?
    normalVelocityNew[:] .+= dropdims(reduce(+, weightsOnEdge[:,:] .* normalVelocity2 .* fEdge2, dims=1), dims=1) # i, iEdge => iEdge
    # TODO boundary condition, will error when any(edgesOnEdge == 0)

    normalVelocity[:] .+= normalVelocityNew[:]

    return nothing
end

function stepSSHParallel!(sshNew, ssh, normalVelocity, bottomDepth, nEdgesOnCell, edgesOnCell, cellsOnCell, areaCell, edgeSignOnCell, dvEdge, dt)

    sshNew[:] .= ssh[:]

    # sum flux through each edge of cell

    bottomDepthOnCellNeighbors = map(iCell -> iCell != 0 ? bottomDepth[iCell] : 0, cellsOnCell)
    normalVelocityOnCellEdges = map(iEdge -> iEdge != 0 ? normalVelocity[iEdge] : 0, edgesOnCell)
    dvEdgeOnCellEdges = map(iEdge -> iEdge != 0 ? dvEdge[iEdge] : 0, edgesOnCell)

    mean_depth = ( bottomDepthOnCellNeighbors .+ reshape(bottomDepth, (1, size(bottomDepth)...) ) ) ./ 2.0
    flux = mean_depth .* normalVelocityOnCellEdges
    sshNew[:] .+= dropdims(reduce(+, flux .* edgeSignOnCell .* dvEdgeOnCellEdges ./ reshape(areaCell, (1, size(areaCell)...)), dims=1), dims=1) # i, iCell => iCell
    # TODO: boundary conditions, will error currently when edgeID == 0

    ssh[:] = sshNew[:]

    return nothing
end

function forwardBackwardStepParallel!(myMPAS_O, dt)

#     cellsOnCell = []
#     for iCell in 1:myMPAS_O.nCells

    stepNormalVelocityParallel!(myMPAS_O.normalVelocityNew,
                                                        myMPAS_O.normalVelocityCurrent,
                                                        myMPAS_O.sshCurrent,
                                                        myMPAS_O.cellsOnEdge,
                                                        myMPAS_O.nEdgesOnEdge,
                                                        myMPAS_O.edgesOnEdge,
                                                        myMPAS_O.weightsOnEdge,
                                                        myMPAS_O.fEdge,
                                                        myMPAS_O.dcEdge,
                                                        dt)

    stepSSHParallel!(myMPAS_O.sshNew,
                                             myMPAS_O.sshCurrent,
                                             myMPAS_O.normalVelocityNew,
                                             myMPAS_O.bottomDepth,
                                             myMPAS_O.nEdgesOnCell,
                                             myMPAS_O.edgesOnCell,
                                             myMPAS_O.cellsOnCell,
                                             myMPAS_O.areaCell,
                                             myMPAS_O.edgeSignOnCell,
                                             myMPAS_O.dvEdge,
                                             dt)
end

function stepNormalVelocity!(normalVelocityNew, normalVelocity, ssh, cellsOnEdge, nEdgesOnEdge, edgesOnEdge, weightsOnEdge, fEdge, dcEdge, dt)
    iEdge = threadIdx().x

    normalVelocityNew[iEdge] = normalVelocity[iEdge]

    # take gradient of sshCurrent across edge
    cell1Index = cellsOnEdge[1,iEdge]
    cell2Index = cellsOnEdge[2,iEdge]

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

    normalVelocity[iEdge] = normalVelocityNew[iEdge]
#     @inbounds myMPAS_O.normalVelocityNew[iEdge] = myMPAS_O.normalVelocityCurrent[iEdge] + dt* normalVelocityTendency#ComputeNormalVelocityTendency(iEdge, myMPAS_O, normalVelocity, ssh)

    return
end

function stepSSH!(sshNew, ssh, normalVelocity, bottomDepth, nEdgesOnCell, edgesOnCell, cellsOnCell, areaCell, edgeSignOnCell, dvEdge, dt)
    iCell = threadIdx().x

    sshNew[iCell] = ssh[iCell]

    # sum flux through each edge of cell
    for i in 1:nEdgesOnCell[iCell]
        edgeID = edgesOnCell[i,iCell]
        neighborCellID = cellsOnCell[i,iCell]

        if neighborCellID !== 0 && edgeID !== 0 # boundary conditions: zero flux at boundaries
            mean_depth = ( bottomDepth[neighborCellID] + bottomDepth[iCell] ) / 2.0
            flux = mean_depth * normalVelocity[edgeID]
            sshNew[iCell] += flux * edgeSignOnCell[iCell,i] * dvEdge[edgeID] / areaCell[iCell]
        end
    end

    ssh[iCell] = sshNew[iCell]

    return
end

function makeStepKernels(myMPAS_O)
    stepNormalVelocityKernel! = stepNormalVelocity!(myMPAS_O.device, myMPAS_O.workGroupSize)
    stepSSHKernel!            = stepSSH!(myMPAS_O.device, myMPAS_O.workGroupSize)
    return stepNormalVelocityKernel!, stepSSHKernel!
end

function forwardBackwardStep!(myMPAS_O, dt, stepNormalVelocityKernel!, stepSSHKernel!)
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

    stepSSHEvent            = stepSSHKernel!(myMPAS_O.sshNew,
                                             myMPAS_O.sshCurrent,
                                             myMPAS_O.normalVelocityNew,
                                             myMPAS_O.bottomDepth,
                                             myMPAS_O.nEdgesOnCell,
                                             myMPAS_O.edgesOnCell,
                                             myMPAS_O.cellsOnCell,
                                             myMPAS_O.areaCell,
                                             myMPAS_O.edgeSignOnCell,
                                             myMPAS_O.dvEdge,
                                             dt, ndrange=myMPAS_O.nCells)
    wait(myMPAS_O.device, stepSSHEvent)
end

function forwardBackwardStep2!(myMPAS_O, dt)
    @cuda threads=myMPAS_O.nEdges stepNormalVelocity!(myMPAS_O.normalVelocityNew,
                                                        myMPAS_O.normalVelocityCurrent,
                                                        myMPAS_O.sshCurrent,
                                                        myMPAS_O.cellsOnEdge,
                                                        myMPAS_O.nEdgesOnEdge,
                                                        myMPAS_O.edgesOnEdge,
                                                        myMPAS_O.weightsOnEdge,
                                                        myMPAS_O.fEdge,
                                                        myMPAS_O.dcEdge,
                                                        dt)

    @cuda threads=myMPAS_O.nCells stepSSH!(myMPAS_O.sshNew,
                                             myMPAS_O.sshCurrent,
                                             myMPAS_O.normalVelocityNew,
                                             myMPAS_O.bottomDepth,
                                             myMPAS_O.nEdgesOnCell,
                                             myMPAS_O.edgesOnCell,
                                             myMPAS_O.cellsOnCell,
                                             myMPAS_O.areaCell,
                                             myMPAS_O.edgeSignOnCell,
                                             myMPAS_O.dvEdge,
                                             dt)
end




myMPAS_O = MPAS_Ocean(false, CODE_ROOT * "MPAS_O_Shallow_Water/Mesh+Initial_Condition+Registry_Files/Periodic",
                "base_mesh.nc", "mesh.nc", periodicity="Periodic"
)
typeof(myMPAS_O)

moveToDevice!(myMPAS_O, GPU())



function gaussianInit!(myMPAS_O)
    sigmax = myMPAS_O.lX/10
    sigmay = myMPAS_O.lY/10

    mux = myMPAS_O.lX/2
    muy = myMPAS_O.lY/2

    for iCell in 1:myMPAS_O.nCells
        myMPAS_O.sshCurrent[iCell] = exp( - (myMPAS_O.xCell[iCell] - mux)^2 / sigmax^2 - (myMPAS_O.yCell[iCell] - muy)^2 / sigmay^2 )
    end

    myMPAS_O.normalVelocityCurrent .= 0

    return nothing
end

gaussianInit!(myMPAS_O)

fig, ax, cbar, col = heatMapMesh(myMPAS_O, myMPAS_O.sshCurrent)
fig

typeof(myMPAS_O.edgeSignOnCell)

moveToDevice!(myMPAS_O, GPU())

# stepNormVelKernel!, stepSSHKernel! = makeStepKernels(myMPAS_O)

nFrames = 30

sshOverTime = zeros(nFrames, myMPAS_O.nCells)

for i in 1:nFrames
    for j in 1:10
        bruh throw an error here!!!
        forwardBackwardStepParallel!(myMPAS_O, 1.0)#, stepNormVelKernel!, stepSSHKernel!)
    end
    sshOverTime[i,:] = myMPAS_O.sshCurrent[:]
end



# stepNormVelKernel!, stepSSHKernel! = makeStepKernels(myMPAS_O)

nFrames = 30

sshOverTime = zeros(nFrames, myMPAS_O.nCells)

for i in 1:nFrames
    for j in 1:10
        forwardBackwardStep2!(myMPAS_O, 1.0)#, stepNormVelKernel!, stepSSHKernel!)
    end
    sshOverTime[i,:] = myMPAS_O.sshCurrent[:]
end

cMax = maximum(abs.(sshOverTime))

fig, ax, _, col = heatMapMesh(myMPAS_O, sshOverTime[1,:], cMin=-cMax, cMax=cMax)

function nextFrame(i)
    col.set_array(sshOverTime[i+1,:])
end

anim = animation.FuncAnimation(fig, nextFrame, frames=nFrames)

ipydisplay.HTML(anim.to_html5_video())

using BenchmarkTools

moveToDevice!(myMPAS_O, CPU())
gaussianInit!(myMPAS_O)
initialSSH[:] = myMPAS_O.sshCurrent[:]
initialNormalVelocity[:] = myMPAS_O.normalVelocityCurrent[:]

function simulate100GPU()
    myMPAS_O.sshCurrent[:] = initialSSH[:]
    myMPAS_O.normalVelocityCurrent[:] = initialNormalVelocity[:]
    moveToDevice!(myMPAS_O, GPU())
    for j in 1:100
        forwardBackwardStep!(myMPAS_O, 1.0, stepNormVelKernel!, stepSSHKernel!)
    end
end

function simulate100CPU()
    gaussianInit!(myMPAS_O)
    moveToDevice!(myMPAS_O, CPU())
    for j in 1:100
        forwardBackward!(myMPAS_O, 1.0, 1.0)
    end
end

@btime simulate100GPU()

@btime simulate100CPU()

SSHTendency = CUDA.zeros(myMPAS_O.nCells)

normalVelocity = CuArray(myMPAS_O.normalVelocityCurrent)
ssh = CuArray(myMPAS_O.sshCurrent)
nEdgesOnCell = CuArray(myMPAS_O.nEdgesOnCell)
edgesOnCell = CuArray(myMPAS_O.edgesOnCell)
cellsOnCell = CuArray(myMPAS_O.cellsOnCell)
bottomDepth = CuArray(myMPAS_O.bottomDepth)
edgeSignOnCell = CuArray(myMPAS_O.edgeSignOnCell)
dvEdge = CuArray(myMPAS_O.dvEdge)
areaCell = CuArray(myMPAS_O.areaCell)

kernel = ComputeSSHTendency(GPU(), 16)

event = kernel(SSHTendency, normalVelocity, ssh, nEdgesOnCell, edgesOnCell,
               cellsOnCell, bottomDepth, edgeSignOnCell, dvEdge, areaCell,
                ndrange=size(SSHTendency))

wait(event)

sshTendency
