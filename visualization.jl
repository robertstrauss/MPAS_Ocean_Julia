## Function for displaying the sea surface height (or other cell-defined variable) on the mesh
# this is not trivial like it is when the grid is a uniform quad, where a simple imshow will work.
using PyCall

mplpatches      = pyimport("matplotlib.patches")
mplcollections  = pyimport("matplotlib.collections")

function heatMapMesh(myMPAS_O, phi; fig="new", ax="new", cbar="new", cmap="seismic", cMin="auto", cMax="auto")
    if fig == "new"
        fig = figure()
    end
    if ax == "new"
        ax = fig.add_subplot(1,1,1)
    end
    if cMin == "auto" && cMax == "auto"
        cMax = maximum(abs.(phi))
        cMin = -cMax
    elseif cMin == "auto"
        cMin = minimum(phi)
    elseif cMax == "auto"
        cMax = maximum(phi)
    end
    m = PyPlot.cm.ScalarMappable(cmap=cmap)
    m.set_clim(cMin, cMax)
    if cbar == "new"
        cbar = colorbar(m, ax=ax, fraction=0.046, pad=0.04)
    end


    xMin = 0.0
    xMax = myMPAS_O.lX + myMPAS_O.gridSpacingMagnitude/2.0
    yMin = 0.0
    yMax = myMPAS_O.lY + myMPAS_O.gridSpacingMagnitude/(2.0*sqrt(3.0))
    ax.axis([xMin,xMax,yMin,yMax])

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

            # edge cases
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
        cellpatch = mplpatches.Polygon(vertices, closed=false)
        append!(patches, [cellpatch])
    end
    localPatches = mplcollections.PatchCollection(patches, cmap=cmap, alpha=1.0)
    localPatches.set_array(phi)
    localPatches.set_clim([cMin, cMax])
    ax.add_collection(localPatches)

    close()
    return fig, ax, cbar, localPatches
end









function edgeHeatMapMesh(myMPAS_O, phi; fig="new", ax="new", cbar="new", cmap="seismic", cMin="auto", cMax="auto", linewidth=1)
    if fig == "new"
        fig = figure()
    end
    if ax == "new"
        ax = fig.add_subplot(1,1,1)
    end
    if cMin == "auto" && cMax == "auto"
        cMax = maximum(abs.(phi))
        cMin = -cMax
    elseif cMin == "auto"
        cMin = minimum(phi)
    elseif cMax == "auto"
        cMax = maximum(phi)
    end
    m = PyPlot.cm.ScalarMappable(cmap=cmap)
    m.set_clim(cMin, cMax)
    if cbar == "new"
        cbar = colorbar(m, ax=ax, fraction=0.046, pad=0.04)
    end

    xMin = 0.0
    xMax = myMPAS_O.lX + myMPAS_O.gridSpacingMagnitude/2.0
    yMin = 0.0
    yMax = myMPAS_O.lY + myMPAS_O.gridSpacingMagnitude/(2.0*sqrt(3.0))
    ax.axis([xMin,xMax,yMin,yMax])

    aspect_ratio = (xMax - xMin)/(yMax - yMin)
    ax.set_aspect(aspect_ratio,adjustable="box")


    patches = []

    for iEdge in 1:myMPAS_O.nEdges

        vertexIndices = myMPAS_O.verticesOnEdge[:,iEdge]

        vertices = [
            myMPAS_O.xVertex[vertexIndices[1]]   myMPAS_O.yVertex[vertexIndices[1]];
            myMPAS_O.xVertex[vertexIndices[2]]   myMPAS_O.yVertex[vertexIndices[2]]
        ]

        # edge cases
        if abs(vertices[1,1] - vertices[2,1]) > 1.01*myMPAS_O.gridSpacingMagnitude
            vertices[1,1] = vertices[2,1] - sign(vertices[1,1] - vertices[2,1]) * myMPAS_O.gridSpacingMagnitude * sqrt(3)/4
        end
        if abs(vertices[1,2] - vertices[2,2]) > 1.01*myMPAS_O.gridSpacingMagnitude
            vertices[2,2] = vertices[1,2] + sign(vertices[1,2] - vertices[2,2]) * myMPAS_O.gridSpacingMagnitude / (2*sqrt(3))
        end

        edgepatch = mplpatches.Polygon(vertices, closed=false, edgecolor=m.to_rgba(phi[iEdge]), linewidth=linewidth)
        append!(patches, [edgepatch])

    end

    localPatches = mplcollections.PatchCollection(patches, cmap=cmap, match_original=true)

    ax.add_collection(localPatches)

    close()
    return fig, ax, cbar, localPatches
end







function signedEdgeHeatMapMesh(myMPAS_O, phi; fig="new", ax="new", cbar="new", cmap="seismic", cMin="auto", cMax="auto", linewidth=1, nudge=0.2)
    if fig == "new"
        fig = figure()
    end
    if ax == "new"
        ax = fig.add_subplot(1,1,1)
    end
    if cMin == "auto" && cMax == "auto"
        cMax = maximum(abs.(phi))
        cMin = -cMax
    elseif cMin == "auto"
        cMin = minimum(phi)
    elseif cMax == "auto"
        cMax = maximum(phi)
    end
    m = PyPlot.cm.ScalarMappable(cmap=cmap)
    m.set_clim(cMin, cMax)
    if cbar == "new"
        cbar = colorbar(m, ax=ax, fraction=0.046, pad=0.04)
    end

    xMin = 0.0
    xMax = myMPAS_O.lX + myMPAS_O.gridSpacingMagnitude/2.0
    yMin = 0.0
    yMax = myMPAS_O.lY + myMPAS_O.gridSpacingMagnitude/(2.0*sqrt(3.0))
    ax.axis([xMin,xMax,yMin,yMax])

    aspect_ratio = (xMax - xMin)/(yMax - yMin)
    ax.set_aspect(aspect_ratio,adjustable="box")


    patches = []

    for iCell in 1:myMPAS_O.nCells
        for i in 1:myMPAS_O.nEdgesOnCell[iCell]
            iEdge = myMPAS_O.edgesOnCell[i, iCell]

            vertexIndices = myMPAS_O.verticesOnEdge[:,iEdge]

            vertices = [
                myMPAS_O.xVertex[vertexIndices[1]]   myMPAS_O.yVertex[vertexIndices[1]];
                myMPAS_O.xVertex[vertexIndices[2]]   myMPAS_O.yVertex[vertexIndices[2]]
            ]


            vertices[:,1] .*= 1.0 - nudge
            vertices[:,1] .+= myMPAS_O.xCell[iCell]*nudge
            vertices[:,2] .*= 1.0 - nudge
            vertices[:,2] .+= myMPAS_O.yCell[iCell]*nudge

            # edge cases
            if abs(vertices[1,1] - vertices[2,1]) > 1.01*myMPAS_O.gridSpacingMagnitude
                vertices[1,:] .= vertices[2,:] #- sign(vertices[1,1] - vertices[2,1]) * myMPAS_O.gridSpacingMagnitude * sqrt(3)/4
            end
            if abs(vertices[1,2] - vertices[2,2]) > 1.01*myMPAS_O.gridSpacingMagnitude
                vertices[2,:] .= vertices[1,:] #+ sign(vertices[1,2] - vertices[2,2]) * myMPAS_O.gridSpacingMagnitude / (2*sqrt(3))
            end

            edgepatch = mplpatches.Polygon(vertices, closed=false, edgecolor=m.to_rgba(phi[iEdge]*myMPAS_O.edgeSignOnCell[iCell, i]), linewidth=linewidth)
            append!(patches, [edgepatch])
        end
    end

    localPatches = mplcollections.PatchCollection(patches, cmap=cmap, match_original=true)

    ax.add_collection(localPatches)

    close()
    return fig, ax, cbar, localPatches
end












function vectorPlotMesh(myMPAS_O, phi; fig="new", ax="new", cbar="new", cmap="seismic", cMin="auto", cMax="auto", linewidth=1)
    if fig == "new"
        fig = figure()
    end
    if ax == "new"
        ax = fig.add_subplot(1,1,1)
    end
    if cMin == "auto" && cMax == "auto"
        cMax = maximum(abs.(phi))
        cMin = -cMax
    elseif cMin == "auto"
        cMin = minimum(phi)
    elseif cMax == "auto"
        cMax = maximum(phi)
    end
    m = PyPlot.cm.ScalarMappable(cmap=cmap)
    m.set_clim(cMin, cMax)
    if cbar == "new"
        cbar = colorbar(m, ax=ax, fraction=0.046, pad=0.04)
    end

    xMin = 0.0
    xMax = myMPAS_O.lX + myMPAS_O.gridSpacingMagnitude/2.0
    yMin = 0.0
    yMax = myMPAS_O.lY + myMPAS_O.gridSpacingMagnitude/(2.0*sqrt(3.0))
    ax.axis([xMin,xMax,yMin,yMax])

    aspect_ratio = (xMax - xMin)/(yMax - yMin)
    ax.set_aspect(aspect_ratio,adjustable="box")

    patches = []

    for iEdge in 1:myMPAS_O.nEdges

        vertexIndices = myMPAS_O.verticesOnEdge[:,iEdge]

        vertices = [
            myMPAS_O.xVertex[vertexIndices[1]]   myMPAS_O.yVertex[vertexIndices[1]];
            myMPAS_O.xVertex[vertexIndices[2]]   myMPAS_O.yVertex[vertexIndices[2]]
        ]

        # edge cases
        if abs(vertices[1,1] - vertices[2,1]) > 1.01*myMPAS_O.gridSpacingMagnitude
            vertices[1,1] = vertices[2,1] - sign(vertices[1,1] - vertices[2,1]) * myMPAS_O.gridSpacingMagnitude * sqrt(3)/4
        end
        if abs(vertices[1,2] - vertices[2,2]) > 1.01*myMPAS_O.gridSpacingMagnitude
            vertices[2,2] = vertices[1,2] + sign(vertices[1,2] - vertices[2,2]) * myMPAS_O.gridSpacingMagnitude / (2*sqrt(3))
        end

        edgepatch = mplpatches.Polygon(vertices, closed=false, edgecolor="k", linewidth=0.5)
        append!(patches, [edgepatch])

    end

    localPatches = mplcollections.PatchCollection(patches, cmap=cmap, match_original=true)

    ax.add_collection(localPatches)

    ax.quiver(myMPAS_O.xEdge, myMPAS_O.yEdge, phi.*cos.(myMPAS_O.angleEdge), phi.*sin.(myMPAS_O.angleEdge))

    close()
    return fig, ax, cbar#, localPatches
end








function vectorStreamPlotMesh(myMPAS_O, phi; fig="new", ax="new", cbar="new", cmap="seismic", cMin="auto", cMax="auto", linecolor="k", linewidth=1.0, density=2.0, xbins=15, ybins=15)
    if fig == "new"
        fig = figure()
    end
    if ax == "new"
        ax = fig.add_subplot(1,1,1)
    end
    if cMin == "auto" && cMax == "auto"
        cMax = maximum(abs.(phi))
        cMin = -cMax
    elseif cMin == "auto"
        cMin = minimum(phi)
    elseif cMax == "auto"
        cMax = maximum(phi)
    end
    m = PyPlot.cm.ScalarMappable(cmap=cmap)
    m.set_clim(cMin, cMax)
    if cbar == "new"
        cbar = colorbar(m, ax=ax, fraction=0.046, pad=0.04)
    end

    xMin = 0.0
    xMax = myMPAS_O.lX + myMPAS_O.gridSpacingMagnitude/2.0
    yMin = 0.0
    yMax = myMPAS_O.lY + myMPAS_O.gridSpacingMagnitude/(2.0*sqrt(3.0))
    ax.axis([xMin,xMax,yMin,yMax])

    aspect_ratio = (xMax - xMin)/(yMax - yMin)
    ax.set_aspect(aspect_ratio,adjustable="box")

#     patches = []

#     for iEdge in 1:myMPAS_O.nEdges

#         vertexIndices = myMPAS_O.verticesOnEdge[:,iEdge]

#         vertices = [
#             myMPAS_O.xVertex[vertexIndices[1]]   myMPAS_O.yVertex[vertexIndices[1]];
#             myMPAS_O.xVertex[vertexIndices[2]]   myMPAS_O.yVertex[vertexIndices[2]]
#         ]

#         # edge cases
#         if abs(vertices[1,1] - vertices[2,1]) > 1.01*myMPAS_O.gridSpacingMagnitude
#             vertices[1,1] = vertices[2,1] - sign(vertices[1,1] - vertices[2,1]) * myMPAS_O.gridSpacingMagnitude * sqrt(3)/4
#         end
#         if abs(vertices[1,2] - vertices[2,2]) > 1.01*myMPAS_O.gridSpacingMagnitude
#             vertices[2,2] = vertices[1,2] + sign(vertices[1,2] - vertices[2,2]) * myMPAS_O.gridSpacingMagnitude / (2*sqrt(3))
#         end

#         edgepatch = mplpatches.Polygon(vertices, closed=false, edgecolor="k", linewidth=0.5)
#         append!(patches, [edgepatch])

#     end

#     localPatches = mplcollections.PatchCollection(patches, cmap=cmap, match_original=true)

#     ax.add_collection(localPatches)


    X = collect(range(0, myMPAS_O.lX, length=xbins))
    Y = collect(range(0, myMPAS_O.lY, length=ybins))

    phiXUniformBinned = zeros(Float64, length(Y), length(X))
    phiYUniformBinned = zeros(Float64, length(Y), length(X))

    for iEdge = 1:myMPAS_O.nEdges
        iBin = argmin(abs.(X .- myMPAS_O.xEdge[iEdge]))
        jBin = argmin(abs.(Y .- myMPAS_O.yEdge[iEdge]))

        phiXUniformBinned[jBin, iBin] += phi[iEdge] * cos(myMPAS_O.angleEdge[iEdge])
        phiYUniformBinned[jBin, iBin] += phi[iEdge] * sin(myMPAS_O.angleEdge[iEdge])
    end

    ax.streamplot(X, Y, phiXUniformBinned, phiYUniformBinned, color=linecolor, density=density, linewidth=linewidth)

    close()

    return fig, ax, cbar#, localPatches
end






function vertexHeatMapMesh(myMPAS_O, phi; fig="new", ax="new", cbar="new", cmap="seismic", cMin="auto", cMax="auto", dotsize=3)
    if fig == "new"
        fig = figure()
    end
    if ax == "new"
        ax = fig.add_subplot(1,1,1)
    end
    if cMin == "auto" && cMax == "auto"
        cMax = maximum(abs.(phi))
        cMin = -cMax
    elseif cMin == "auto"
        cMin = minimum(phi)
    elseif cMax == "auto"
        cMax = maximum(phi)
    end
    m = PyPlot.cm.ScalarMappable(cmap=cmap)
    m.set_clim(cMin, cMax)
    if cbar == "new"
        cbar = colorbar(m, ax=ax, fraction=0.046, pad=0.04)
    end

    xMin = 0.0
    xMax = myMPAS_O.lX + myMPAS_O.gridSpacingMagnitude/2.0
    yMin = 0.0
    yMax = myMPAS_O.lY + myMPAS_O.gridSpacingMagnitude/(2.0*sqrt(3.0))
    ax.axis([xMin,xMax,yMin,yMax])

    aspect_ratio = (xMax - xMin)/(yMax - yMin)
    ax.set_aspect(aspect_ratio,adjustable="box")


    patches = []

    for iEdge in 1:myMPAS_O.nEdges

        vertexIndices = myMPAS_O.verticesOnEdge[:,iEdge]

        vertices = [
            myMPAS_O.xVertex[vertexIndices[1]]   myMPAS_O.yVertex[vertexIndices[1]];
            myMPAS_O.xVertex[vertexIndices[2]]   myMPAS_O.yVertex[vertexIndices[2]]
        ]

        # edge cases
        if abs(vertices[1,1] - vertices[2,1]) > 1.01*myMPAS_O.gridSpacingMagnitude
            vertices[1,1] = vertices[2,1] - sign(vertices[1,1] - vertices[2,1]) * myMPAS_O.gridSpacingMagnitude * sqrt(3)/4
        end
        if abs(vertices[1,2] - vertices[2,2]) > 1.01*myMPAS_O.gridSpacingMagnitude
            vertices[2,2] = vertices[1,2] + sign(vertices[1,2] - vertices[2,2]) * myMPAS_O.gridSpacingMagnitude / (2*sqrt(3))
        end

        edgepatch = mplpatches.Polygon(vertices, closed=false, edgecolor="k", linewidth=0.5)
        append!(patches, [edgepatch])

    end

    localPatches = mplcollections.PatchCollection(patches, cmap=cmap, match_original=true)

    ax.add_collection(localPatches)

    ax.scatter(myMPAS_O.xVertex, myMPAS_O.yVertex, s=dotsize, c=phi, cmap=cmap)

    close()
    return fig, ax, cbar, localPatches
end
