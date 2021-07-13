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
#     if ( typeof(getfield(mpasOcean, field)) <: CuArray )
#         append!(fields, [field])
#         append!(values, [getfield(mpasOcean, field)])
#     end
# #     @show typeof(getfield(mpasOcean, field))
# end
#
#
# append!(fields, [:device])
# append!(values, [mpasOcean.device])
#
# append!(fields, [:workGroupSize])
# append!(values, [mpasOcean.workGroupSize])
#
# append!(fields, [:nEdges])
# append!(values, [mpasOcean.nEdges])
#
# append!(fields, [:nCells])
# append!(values, [mpasOcean.nCells])
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

function forwardBackwardStepParallel!(mpasOcean, dt)

#     cellsOnCell = []
#     for iCell in 1:mpasOcean.nCells

    stepNormalVelocityParallel!(mpasOcean.normalVelocityNew,
                                                        mpasOcean.normalVelocityCurrent,
                                                        mpasOcean.sshCurrent,
                                                        mpasOcean.cellsOnEdge,
                                                        mpasOcean.nEdgesOnEdge,
                                                        mpasOcean.edgesOnEdge,
                                                        mpasOcean.weightsOnEdge,
                                                        mpasOcean.fEdge,
                                                        mpasOcean.dcEdge,
                                                        dt)

    stepSSHParallel!(mpasOcean.sshNew,
                                             mpasOcean.sshCurrent,
                                             mpasOcean.normalVelocityNew,
                                             mpasOcean.bottomDepth,
                                             mpasOcean.nEdgesOnCell,
                                             mpasOcean.edgesOnCell,
                                             mpasOcean.cellsOnCell,
                                             mpasOcean.areaCell,
                                             mpasOcean.edgeSignOnCell,
                                             mpasOcean.dvEdge,
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
#     @inbounds mpasOcean.normalVelocityNew[iEdge] = mpasOcean.normalVelocityCurrent[iEdge] + dt* normalVelocityTendency#ComputeNormalVelocityTendency(iEdge, mpasOcean, normalVelocity, ssh)

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

function makeStepKernels(mpasOcean)
    stepNormalVelocityKernel! = stepNormalVelocity!(mpasOcean.device, mpasOcean.workGroupSize)
    stepSSHKernel!            = stepSSH!(mpasOcean.device, mpasOcean.workGroupSize)
    return stepNormalVelocityKernel!, stepSSHKernel!
end

function forwardBackwardStep!(mpasOcean, dt, stepNormalVelocityKernel!, stepSSHKernel!)
    stepNormalVelocityEvent = stepNormalVelocityKernel!(mpasOcean.normalVelocityNew,
                                                        mpasOcean.normalVelocityCurrent,
                                                        mpasOcean.sshCurrent,
                                                        mpasOcean.cellsOnEdge,
                                                        mpasOcean.nEdgesOnEdge,
                                                        mpasOcean.edgesOnEdge,
                                                        mpasOcean.weightsOnEdge,
                                                        mpasOcean.fEdge,
                                                        mpasOcean.dcEdge,
                                                        dt, ndrange=mpasOcean.nEdges)
    wait(mpasOcean.device, stepNormalVelocityEvent)

    stepSSHEvent            = stepSSHKernel!(mpasOcean.sshNew,
                                             mpasOcean.sshCurrent,
                                             mpasOcean.normalVelocityNew,
                                             mpasOcean.bottomDepth,
                                             mpasOcean.nEdgesOnCell,
                                             mpasOcean.edgesOnCell,
                                             mpasOcean.cellsOnCell,
                                             mpasOcean.areaCell,
                                             mpasOcean.edgeSignOnCell,
                                             mpasOcean.dvEdge,
                                             dt, ndrange=mpasOcean.nCells)
    wait(mpasOcean.device, stepSSHEvent)
end

function forwardBackwardStep2!(mpasOcean, dt)
    @cuda threads=mpasOcean.nEdges stepNormalVelocity!(mpasOcean.normalVelocityNew,
                                                        mpasOcean.normalVelocityCurrent,
                                                        mpasOcean.sshCurrent,
                                                        mpasOcean.cellsOnEdge,
                                                        mpasOcean.nEdgesOnEdge,
                                                        mpasOcean.edgesOnEdge,
                                                        mpasOcean.weightsOnEdge,
                                                        mpasOcean.fEdge,
                                                        mpasOcean.dcEdge,
                                                        dt)

    @cuda threads=mpasOcean.nCells stepSSH!(mpasOcean.sshNew,
                                             mpasOcean.sshCurrent,
                                             mpasOcean.normalVelocityNew,
                                             mpasOcean.bottomDepth,
                                             mpasOcean.nEdgesOnCell,
                                             mpasOcean.edgesOnCell,
                                             mpasOcean.cellsOnCell,
                                             mpasOcean.areaCell,
                                             mpasOcean.edgeSignOnCell,
                                             mpasOcean.dvEdge,
                                             dt)
end




mpasOcean = MPAS_Ocean(false, CODE_ROOT * "MPAS_O_Shallow_Water/Mesh+Initial_Condition+Registry_Files/Periodic",
                "base_mesh.nc", "mesh.nc", periodicity="Periodic"
)
typeof(mpasOcean)

moveToDevice!(mpasOcean, GPU())



function gaussianInit!(mpasOcean)
    sigmax = mpasOcean.lX/10
    sigmay = mpasOcean.lY/10

    mux = mpasOcean.lX/2
    muy = mpasOcean.lY/2

    for iCell in 1:mpasOcean.nCells
        mpasOcean.sshCurrent[iCell] = exp( - (mpasOcean.xCell[iCell] - mux)^2 / sigmax^2 - (mpasOcean.yCell[iCell] - muy)^2 / sigmay^2 )
    end

    mpasOcean.normalVelocityCurrent .= 0

    return nothing
end

gaussianInit!(mpasOcean)

fig, ax, cbar, col = heatMapMesh(mpasOcean, mpasOcean.sshCurrent)
fig

typeof(mpasOcean.edgeSignOnCell)

moveToDevice!(mpasOcean, GPU())

# stepNormVelKernel!, stepSSHKernel! = makeStepKernels(mpasOcean)

nFrames = 30

sshOverTime = zeros(nFrames, mpasOcean.nCells)

for i in 1:nFrames
    for j in 1:10
        bruh throw an error here!!!
        forwardBackwardStepParallel!(mpasOcean, 1.0)#, stepNormVelKernel!, stepSSHKernel!)
    end
    sshOverTime[i,:] = mpasOcean.sshCurrent[:]
end



# stepNormVelKernel!, stepSSHKernel! = makeStepKernels(mpasOcean)

nFrames = 30

sshOverTime = zeros(nFrames, mpasOcean.nCells)

for i in 1:nFrames
    for j in 1:10
        forwardBackwardStep2!(mpasOcean, 1.0)#, stepNormVelKernel!, stepSSHKernel!)
    end
    sshOverTime[i,:] = mpasOcean.sshCurrent[:]
end

cMax = maximum(abs.(sshOverTime))

fig, ax, _, col = heatMapMesh(mpasOcean, sshOverTime[1,:], cMin=-cMax, cMax=cMax)

function nextFrame(i)
    col.set_array(sshOverTime[i+1,:])
end

anim = animation.FuncAnimation(fig, nextFrame, frames=nFrames)

ipydisplay.HTML(anim.to_html5_video())

using BenchmarkTools

moveToDevice!(mpasOcean, CPU())
gaussianInit!(mpasOcean)
initialSSH[:] = mpasOcean.sshCurrent[:]
initialNormalVelocity[:] = mpasOcean.normalVelocityCurrent[:]

function simulate100GPU()
    mpasOcean.sshCurrent[:] = initialSSH[:]
    mpasOcean.normalVelocityCurrent[:] = initialNormalVelocity[:]
    moveToDevice!(mpasOcean, GPU())
    for j in 1:100
        forwardBackwardStep!(mpasOcean, 1.0, stepNormVelKernel!, stepSSHKernel!)
    end
end

function simulate100CPU()
    gaussianInit!(mpasOcean)
    moveToDevice!(mpasOcean, CPU())
    for j in 1:100
        forwardBackward!(mpasOcean, 1.0, 1.0)
    end
end

@btime simulate100GPU()

@btime simulate100CPU()

SSHTendency = CUDA.zeros(mpasOcean.nCells)

normalVelocity = CuArray(mpasOcean.normalVelocityCurrent)
ssh = CuArray(mpasOcean.sshCurrent)
nEdgesOnCell = CuArray(mpasOcean.nEdgesOnCell)
edgesOnCell = CuArray(mpasOcean.edgesOnCell)
cellsOnCell = CuArray(mpasOcean.cellsOnCell)
bottomDepth = CuArray(mpasOcean.bottomDepth)
edgeSignOnCell = CuArray(mpasOcean.edgeSignOnCell)
dvEdge = CuArray(mpasOcean.dvEdge)
areaCell = CuArray(mpasOcean.areaCell)

kernel = ComputeSSHTendency(GPU(), 16)

event = kernel(SSHTendency, normalVelocity, ssh, nEdgesOnCell, edgesOnCell,
               cellsOnCell, bottomDepth, edgeSignOnCell, dvEdge, areaCell,
                ndrange=size(SSHTendency))

wait(event)

sshTendency
