include("mode_forward.jl")
include("mode_init.jl")
include("visualization.jl")

using PyPlot
using PyCall

patch = pyimport("matplotlib.patches")
col = pyimport("matplotlib.collections")
animation = pyimport("matplotlib.animation")

testCase1 = true
if testCase1
    mpasOcean = MPAS_Ocean(mesh_directory = "MPAS_O_Shallow_Water/Mesh+Initial_Condition+Registry_Files/Periodic",
                        base_mesh_file_name = "base_mesh.nc",
                        mesh_file_name = "mesh.nc")


    function gaussianInit!(mpasOcean)
        # define gaussian initial condition
        for iCell in 1:mpasOcean.nCells
            mpasOcean.sshCurrent[iCell] = exp(
                (
                    -(mpasOcean.yCell[iCell] - 20000)^2
                    -(mpasOcean.xCell[iCell] - 30000)^2
                ) / ( 1.0e7 )
            )
        end

        mpasOcean.normalVelocityCurrent = zeros(mpasOcean.nEdges)
    end

    # calculate dt based on CFL condition
    dt = 0.1 * minimum(mpasOcean.dcEdge) / sqrt( gravity  * maximum(mpasOcean.bottomDepth) )
    println("dt: ", dt)

    simulationTime = dt*1000

    nFrames = 100

    sshOverTimeFE = zeros(Float64, (nFrames, mpasOcean.nCells))
    sshOverTimeFB = zeros(Float64, (nFrames, mpasOcean.nCells))


    # plot sea surface height every interval of 1.0e-3 seconds
    function cb!(mpasOcean, t, arr)
        if t%(simulationTime/nFrames) <= dt
            frame = 1+Int(floor(t*nFrames/simulationTime))
            if frame <= size(arr)[1]
                arr[frame,:] = mpasOcean.sshCurrent[:]
            end
        end
    end
    cbFE(mpasOcean, t) = cb!(mpasOcean, t, sshOverTimeFE)
    cbFB(mpasOcean, t) = cb!(mpasOcean, t, sshOverTimeFB)

    # solve it, integrate!
    gaussianInit!(mpasOcean)
    forwardEuler!(mpasOcean, dt, simulationTime, cbFE)
    gaussianInit!(mpasOcean)
    forwardBackward!(mpasOcean, dt, simulationTime, cbFB)




    fig = figure()

    ax1 = fig.add_subplot(1, 3, 1)
    ax2 = fig.add_subplot(1, 3, 2)
    ax3 = fig.add_subplot(1, 3, 3)

    cMax1 = maximum(abs.(sshOverTimeFE))
    cMax2 = maximum(abs.(sshOverTimeFB))
    cMax3 = maximum(abs.(sshOverTimeFE .- sshOverTimeFB))

    ccbar1 = "new"
    ccbar2 = "new"
    ccbar3 = "new"

    function nextFrame(i)

        global ccbar1, ccbar2, ccbar3

        j = i+1

        ffig1, aax1, ccbar1 = heatMapMesh(mpasOcean, sshOverTimeFE[j,:]; fig=fig, ax=ax1, cbar=ccbar1, cMin=-cMax1, cMax=cMax1)
        aax1.set_title("Forward Euler")
        ffig2, aax2, ccbar2 = heatMapMesh(mpasOcean, sshOverTimeFB[j,:]; fig=fig, ax=ax2, cbar=ccbar2, cMin=-cMax2, cMax=cMax2)
        aax2.set_title("Forward Backward")
        ffig3, aax3, ccbar3 = heatMapMesh(mpasOcean, sshOverTimeFE[j,:] .- sshOverTimeFB[j,:]; fig=fig, ax=ax3, cbar=ccbar3, cMin=-cMax3, cMax=cMax3)
        aax3.set_title("Difference")

        fig.suptitle("frame $j of $nFrames")

        return [aax1, aax2, aax3]

    end

    for i in 0:nFrames-1
        nextFrame(i)
        display(fig)
        # sleep(0.01)
    end

    # anim = animation.FuncAnimation(fig, nextFrame)
    #
    # anim.save("plots/test2/eulervsFB.mp4")

end
