include("mode_forward.jl")
include("mode_init.jl")

using PyPlot
using PyCall

patch = pyimport("matplotlib.patches")
col = pyimport("matplotlib.collections")
animation = pyimport("matplotlib.animation")

function heatMapMesh(myMPAS_O, phi; fig="new", ax="new", cbar="new", cmap="seismic", cMin=-1.0, cMax=1.0)
    if fig == "new"
        fig = figure()
    end
    if ax == "new"
        ax = fig.add_subplot(1,1,1)
    end

    xMin = 0.0
    xMax = myMPAS_O.lX + myMPAS_O.gridSpacingMagnitude/2.0
    yMin = 0.0
    yMax = myMPAS_O.lY + myMPAS_O.gridSpacingMagnitude/(2.0*sqrt(3.0))
    aspect_ratio = (xMax - xMin)/(yMax - yMin)
    ax.set_aspect(aspect_ratio,adjustable="box")

    patches = []
    ComparisonTolerance = 10.0^(-10.0)
    for iCell in 1:myMPAS_O.nCells
        xCell = myMPAS_O.xCell[iCell]
        yCell = myMPAS_O.yCell[iCell]
        nVertices = myMPAS_O.nEdgesOnCell[iCell]
        vertexIndices = myMPAS_O.verticesOnCell[:,iCell]
        vertices = zeros(Float64, (nVertices, 2))
        for iVertexOnCell in 1:nVertices
            xVertex = myMPAS_O.xVertex[vertexIndices[iVertexOnCell]]
            yVertex = myMPAS_O.yVertex[vertexIndices[iVertexOnCell]]
            if abs(yVertex - yCell) > (2.0/sqrt(3.0))*myMPAS_O.gridSpacingMagnitude && yVertex < yCell
                yVertex = yCell + myMPAS_O.gridSpacingMagnitude/sqrt(3.0)
            end
            if abs(yVertex - yCell) > (2.0/sqrt(3.0))*myMPAS_O.gridSpacingMagnitude && yVertex > yCell
                yVertex = yCell - myMPAS_O.gridSpacingMagnitude/sqrt(3.0)
            end
            if abs(xVertex - xCell) > myMPAS_O.gridSpacingMagnitude && xVertex < xCell
                if abs(yVertex - (yCell + myMPAS_O.gridSpacingMagnitude/sqrt(3.0))) < ComparisonTolerance
                    xVertex = xCell
                elseif abs(yVertex - (yCell - myMPAS_O.gridSpacingMagnitude/sqrt(3.0))) < ComparisonTolerance
                    xVertex = xCell
                else
                    xVertex = xCell + 0.5*myMPAS_O.gridSpacingMagnitude
                end
            end
            vertices[iVertexOnCell,1] = xVertex
            vertices[iVertexOnCell,2] = yVertex
        end
        cellpatch = patch.Polygon(vertices, closed=true, alpha=0.5, linestyle="--", linewidth=3, edgecolor="k")
        append!(patches, [cellpatch])
    end
    localPatches = col.PatchCollection(patches, cmap=cmap, alpha=1.0)
    localPatches.set_array(phi)
    localPatches.set_clim([cMin, cMax])
    ax.add_collection(localPatches)
    ax.axis([xMin,xMax,yMin,yMax])

    if cbar == "new"
        m = PyPlot.cm.ScalarMappable(cmap=cmap)
        m.set_clim(cMin, cMax)
        cbar = colorbar(m, ax=ax)
    end

    close()
    return fig, ax, cbar
end


testCase1 = true
if testCase1
    myMPAS_O = MPAS_Ocean(mesh_directory = "MPAS_O_Shallow_Water/Mesh+Initial_Condition+Registry_Files/Periodic",
                        base_mesh_file_name = "base_mesh.nc",
                        mesh_file_name = "mesh.nc")


    function gaussianInit!(myMPAS_O)
        # define gaussian initial condition
        for iCell in 1:myMPAS_O.nCells
            myMPAS_O.sshCurrent[iCell] = exp(
                (
                    -(myMPAS_O.yCell[iCell] - 20000)^2
                    -(myMPAS_O.xCell[iCell] - 30000)^2
                ) / ( 1.0e7 )
            )
        end

        myMPAS_O.normalVelocityCurrent = zeros(myMPAS_O.nEdges)
    end

    # calculate dt based on CFL condition
    dt = 0.1 * minimum(myMPAS_O.dcEdge) / sqrt( gravity  * maximum(myMPAS_O.bottomDepth) )
    println("dt: ", dt)

    simulationTime = dt*1000

    nFrames = 100

    sshOverTimeFE = zeros(Float64, (nFrames, myMPAS_O.nCells))
    sshOverTimeFB = zeros(Float64, (nFrames, myMPAS_O.nCells))


    # plot sea surface height every interval of 1.0e-3 seconds
    function cb!(myMPAS_O, t, arr)
        if t%(simulationTime/nFrames) <= dt
            frame = 1+Int(floor(t*nFrames/simulationTime))
            if frame <= size(arr)[1]
                arr[frame,:] = myMPAS_O.sshCurrent[:]
            end
        end
    end
    cbFE(myMPAS_O, t) = cb!(myMPAS_O, t, sshOverTimeFE)
    cbFB(myMPAS_O, t) = cb!(myMPAS_O, t, sshOverTimeFB)

    # solve it, integrate!
    gaussianInit!(myMPAS_O)
    forwardEuler!(myMPAS_O, dt, simulationTime, cbFE)
    gaussianInit!(myMPAS_O)
    forwardBackward!(myMPAS_O, dt, simulationTime, cbFB)




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

        ffig1, aax1, ccbar1 = heatMapMesh(myMPAS_O, sshOverTimeFE[j,:]; fig=fig, ax=ax1, cbar=ccbar1, cMin=-cMax1, cMax=cMax1)
        aax1.set_title("Forward Euler")
        ffig2, aax2, ccbar2 = heatMapMesh(myMPAS_O, sshOverTimeFB[j,:]; fig=fig, ax=ax2, cbar=ccbar2, cMin=-cMax2, cMax=cMax2)
        aax2.set_title("Forward Backward")
        ffig3, aax3, ccbar3 = heatMapMesh(myMPAS_O, sshOverTimeFE[j,:] .- sshOverTimeFB[j,:]; fig=fig, ax=ax3, cbar=ccbar3, cMin=-cMax3, cMax=cMax3)
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
