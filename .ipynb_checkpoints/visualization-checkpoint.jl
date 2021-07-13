## Function for displaying the sea surface height (or other cell-defined variable) on the mesh
# this is not trivial like it is when the grid is a uniform quad, where a simple imshow will work.
using PyCall
using PyPlot

mplpatches      = pyimport("matplotlib.patches")
mplcollections  = pyimport("matplotlib.collections")

function heatMapMesh(mpasOcean, phi; fig="new", ax="new", cbar="new", cmap="seismic", cMin="auto", cMax="auto")
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
    xMax = mpasOcean.lX + mpasOcean.gridSpacingMagnitude/2.0
    yMin = 0.0
    yMax = mpasOcean.lY + mpasOcean.gridSpacingMagnitude/(2.0*sqrt(3.0))
    ax.axis([xMin,xMax,yMin,yMax])

    aspect_ratio = (xMax - xMin)/(yMax - yMin)
    ax.set_aspect(aspect_ratio,adjustable="box")

    patches = []
    ComparisonTolerance = 10.0^(-10.0)
    for iCell in 1:mpasOcean.nCells
        xCell = mpasOcean.xCell[iCell]
        yCell = mpasOcean.yCell[iCell]
        nVertices = mpasOcean.nEdgesOnCell[iCell]
        vertexIndices = mpasOcean.verticesOnCell[:,iCell]
        vertices = zeros(Float64, (nVertices, 2))
        for iVertexOnCell in 1:nVertices
            xVertex = mpasOcean.xVertex[vertexIndices[iVertexOnCell]]
            yVertex = mpasOcean.yVertex[vertexIndices[iVertexOnCell]]

            # edge cases
            if abs(yVertex - yCell) > (2.0/sqrt(3.0))*mpasOcean.gridSpacingMagnitude && yVertex < yCell
                yVertex = yCell + mpasOcean.gridSpacingMagnitude/sqrt(3.0)
            end
            if abs(yVertex - yCell) > (2.0/sqrt(3.0))*mpasOcean.gridSpacingMagnitude && yVertex > yCell
                yVertex = yCell - mpasOcean.gridSpacingMagnitude/sqrt(3.0)
            end
            if abs(xVertex - xCell) > mpasOcean.gridSpacingMagnitude && xVertex < xCell
                if abs(yVertex - (yCell + mpasOcean.gridSpacingMagnitude/sqrt(3.0))) < ComparisonTolerance
                    xVertex = xCell
                elseif abs(yVertex - (yCell - mpasOcean.gridSpacingMagnitude/sqrt(3.0))) < ComparisonTolerance
                    xVertex = xCell
                else
                    xVertex = xCell + 0.5*mpasOcean.gridSpacingMagnitude
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









function edgeHeatMapMesh(mpasOcean, phi; fig="new", ax="new", cbar="new", cmap="seismic", cMin="auto", cMax="auto", linewidth=1)
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
    xMax = mpasOcean.lX + mpasOcean.gridSpacingMagnitude/2.0
    yMin = 0.0
    yMax = mpasOcean.lY + mpasOcean.gridSpacingMagnitude/(2.0*sqrt(3.0))
    ax.axis([xMin,xMax,yMin,yMax])

    aspect_ratio = (xMax - xMin)/(yMax - yMin)
    ax.set_aspect(aspect_ratio,adjustable="box")


    patches = []

    for iEdge in 1:mpasOcean.nEdges

        vertexIndices = mpasOcean.verticesOnEdge[:,iEdge]

        vertices = [
            mpasOcean.xVertex[vertexIndices[1]]   mpasOcean.yVertex[vertexIndices[1]];
            mpasOcean.xVertex[vertexIndices[2]]   mpasOcean.yVertex[vertexIndices[2]]
        ]

        # edge cases
        if abs(vertices[1,1] - vertices[2,1]) > 1.01*mpasOcean.gridSpacingMagnitude
            vertices[1,1] = vertices[2,1] - sign(vertices[1,1] - vertices[2,1]) * mpasOcean.gridSpacingMagnitude * sqrt(3)/4
        end
        if abs(vertices[1,2] - vertices[2,2]) > 1.01*mpasOcean.gridSpacingMagnitude
            vertices[2,2] = vertices[1,2] + sign(vertices[1,2] - vertices[2,2]) * mpasOcean.gridSpacingMagnitude / (2*sqrt(3))
        end

        edgepatch = mplpatches.Polygon(vertices, closed=false, edgecolor=m.to_rgba(phi[iEdge]), linewidth=linewidth)
        append!(patches, [edgepatch])

    end

    localPatches = mplcollections.PatchCollection(patches, cmap=cmap, match_original=true)

    ax.add_collection(localPatches)

    close()
    return fig, ax, cbar, localPatches
end







function signedEdgeHeatMapMesh(mpasOcean, phi; fig="new", ax="new", cbar="new", cmap="seismic", cMin="auto", cMax="auto", linewidth=1, nudge=0.2)
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
    xMax = mpasOcean.lX + mpasOcean.gridSpacingMagnitude/2.0
    yMin = 0.0
    yMax = mpasOcean.lY + mpasOcean.gridSpacingMagnitude/(2.0*sqrt(3.0))
    ax.axis([xMin,xMax,yMin,yMax])

    aspect_ratio = (xMax - xMin)/(yMax - yMin)
    ax.set_aspect(aspect_ratio,adjustable="box")


    patches = []

    for iCell in 1:mpasOcean.nCells
        for i in 1:mpasOcean.nEdgesOnCell[iCell]
            iEdge = mpasOcean.edgesOnCell[i, iCell]

            vertexIndices = mpasOcean.verticesOnEdge[:,iEdge]

            vertices = [
                mpasOcean.xVertex[vertexIndices[1]]   mpasOcean.yVertex[vertexIndices[1]];
                mpasOcean.xVertex[vertexIndices[2]]   mpasOcean.yVertex[vertexIndices[2]]
            ]


            vertices[:,1] .*= 1.0 - nudge
            vertices[:,1] .+= mpasOcean.xCell[iCell]*nudge
            vertices[:,2] .*= 1.0 - nudge
            vertices[:,2] .+= mpasOcean.yCell[iCell]*nudge

            # edge cases
            if abs(vertices[1,1] - vertices[2,1]) > 1.01*mpasOcean.gridSpacingMagnitude
                vertices[1,:] .= vertices[2,:] #- sign(vertices[1,1] - vertices[2,1]) * mpasOcean.gridSpacingMagnitude * sqrt(3)/4
            end
            if abs(vertices[1,2] - vertices[2,2]) > 1.01*mpasOcean.gridSpacingMagnitude
                vertices[2,:] .= vertices[1,:] #+ sign(vertices[1,2] - vertices[2,2]) * mpasOcean.gridSpacingMagnitude / (2*sqrt(3))
            end

            edgepatch = mplpatches.Polygon(vertices, closed=false, edgecolor=m.to_rgba(phi[iEdge]*mpasOcean.edgeSignOnCell[iCell, i]), linewidth=linewidth)
            append!(patches, [edgepatch])
        end
    end

    localPatches = mplcollections.PatchCollection(patches, cmap=cmap, match_original=true)

    ax.add_collection(localPatches)

    close()
    return fig, ax, cbar, localPatches
end












function vectorPlotMesh(mpasOcean, phi; fig="new", ax="new", cbar="new", cmap="seismic", cMin="auto", cMax="auto", linewidth=1)
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
    xMax = mpasOcean.lX + mpasOcean.gridSpacingMagnitude/2.0
    yMin = 0.0
    yMax = mpasOcean.lY + mpasOcean.gridSpacingMagnitude/(2.0*sqrt(3.0))
    ax.axis([xMin,xMax,yMin,yMax])

    aspect_ratio = (xMax - xMin)/(yMax - yMin)
    ax.set_aspect(aspect_ratio,adjustable="box")

    patches = []

    for iEdge in 1:mpasOcean.nEdges

        vertexIndices = mpasOcean.verticesOnEdge[:,iEdge]

        vertices = [
            mpasOcean.xVertex[vertexIndices[1]]   mpasOcean.yVertex[vertexIndices[1]];
            mpasOcean.xVertex[vertexIndices[2]]   mpasOcean.yVertex[vertexIndices[2]]
        ]

        # edge cases
        if abs(vertices[1,1] - vertices[2,1]) > 1.01*mpasOcean.gridSpacingMagnitude
            vertices[1,1] = vertices[2,1] - sign(vertices[1,1] - vertices[2,1]) * mpasOcean.gridSpacingMagnitude * sqrt(3)/4
        end
        if abs(vertices[1,2] - vertices[2,2]) > 1.01*mpasOcean.gridSpacingMagnitude
            vertices[2,2] = vertices[1,2] + sign(vertices[1,2] - vertices[2,2]) * mpasOcean.gridSpacingMagnitude / (2*sqrt(3))
        end

        edgepatch = mplpatches.Polygon(vertices, closed=false, edgecolor="k", linewidth=0.5)
        append!(patches, [edgepatch])

    end

    localPatches = mplcollections.PatchCollection(patches, cmap=cmap, match_original=true)

    ax.add_collection(localPatches)

    ax.quiver(mpasOcean.xEdge, mpasOcean.yEdge, phi.*cos.(mpasOcean.angleEdge), phi.*sin.(mpasOcean.angleEdge))

    close()
    return fig, ax, cbar#, localPatches
end








function vectorStreamPlotMesh(mpasOcean, phi; fig="new", ax="new", cbar="new", cmap="seismic", cMin="auto", cMax="auto", linecolor="k", linewidth=1.0, density=2.0, xbins=15, ybins=15)
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
    xMax = mpasOcean.lX + mpasOcean.gridSpacingMagnitude/2.0
    yMin = 0.0
    yMax = mpasOcean.lY + mpasOcean.gridSpacingMagnitude/(2.0*sqrt(3.0))
    ax.axis([xMin,xMax,yMin,yMax])

    aspect_ratio = (xMax - xMin)/(yMax - yMin)
    ax.set_aspect(aspect_ratio,adjustable="box")

#     patches = []

#     for iEdge in 1:mpasOcean.nEdges

#         vertexIndices = mpasOcean.verticesOnEdge[:,iEdge]

#         vertices = [
#             mpasOcean.xVertex[vertexIndices[1]]   mpasOcean.yVertex[vertexIndices[1]];
#             mpasOcean.xVertex[vertexIndices[2]]   mpasOcean.yVertex[vertexIndices[2]]
#         ]

#         # edge cases
#         if abs(vertices[1,1] - vertices[2,1]) > 1.01*mpasOcean.gridSpacingMagnitude
#             vertices[1,1] = vertices[2,1] - sign(vertices[1,1] - vertices[2,1]) * mpasOcean.gridSpacingMagnitude * sqrt(3)/4
#         end
#         if abs(vertices[1,2] - vertices[2,2]) > 1.01*mpasOcean.gridSpacingMagnitude
#             vertices[2,2] = vertices[1,2] + sign(vertices[1,2] - vertices[2,2]) * mpasOcean.gridSpacingMagnitude / (2*sqrt(3))
#         end

#         edgepatch = mplpatches.Polygon(vertices, closed=false, edgecolor="k", linewidth=0.5)
#         append!(patches, [edgepatch])

#     end

#     localPatches = mplcollections.PatchCollection(patches, cmap=cmap, match_original=true)

#     ax.add_collection(localPatches)


    X = collect(range(0, mpasOcean.lX, length=xbins))
    Y = collect(range(0, mpasOcean.lY, length=ybins))

    phiXUniformBinned = zeros(Float64, length(Y), length(X))
    phiYUniformBinned = zeros(Float64, length(Y), length(X))

    for iEdge = 1:mpasOcean.nEdges
        iBin = argmin(abs.(X .- mpasOcean.xEdge[iEdge]))
        jBin = argmin(abs.(Y .- mpasOcean.yEdge[iEdge]))

        phiXUniformBinned[jBin, iBin] += phi[iEdge] * cos(mpasOcean.angleEdge[iEdge])
        phiYUniformBinned[jBin, iBin] += phi[iEdge] * sin(mpasOcean.angleEdge[iEdge])
    end

    ax.streamplot(X, Y, phiXUniformBinned, phiYUniformBinned, color=linecolor, density=density, linewidth=linewidth)

    close()

    return fig, ax, cbar#, localPatches
end






function vertexHeatMapMesh(mpasOcean, phi; fig="new", ax="new", cbar="new", cmap="seismic", cMin="auto", cMax="auto", dotsize=3)
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
    xMax = mpasOcean.lX + mpasOcean.gridSpacingMagnitude/2.0
    yMin = 0.0
    yMax = mpasOcean.lY + mpasOcean.gridSpacingMagnitude/(2.0*sqrt(3.0))
    ax.axis([xMin,xMax,yMin,yMax])

    aspect_ratio = (xMax - xMin)/(yMax - yMin)
    ax.set_aspect(aspect_ratio,adjustable="box")


    patches = []

    for iEdge in 1:mpasOcean.nEdges

        vertexIndices = mpasOcean.verticesOnEdge[:,iEdge]

        vertices = [
            mpasOcean.xVertex[vertexIndices[1]]   mpasOcean.yVertex[vertexIndices[1]];
            mpasOcean.xVertex[vertexIndices[2]]   mpasOcean.yVertex[vertexIndices[2]]
        ]

        # edge cases
        if abs(vertices[1,1] - vertices[2,1]) > 1.01*mpasOcean.gridSpacingMagnitude
            vertices[1,1] = vertices[2,1] - sign(vertices[1,1] - vertices[2,1]) * mpasOcean.gridSpacingMagnitude * sqrt(3)/4
        end
        if abs(vertices[1,2] - vertices[2,2]) > 1.01*mpasOcean.gridSpacingMagnitude
            vertices[2,2] = vertices[1,2] + sign(vertices[1,2] - vertices[2,2]) * mpasOcean.gridSpacingMagnitude / (2*sqrt(3))
        end

        edgepatch = mplpatches.Polygon(vertices, closed=false, edgecolor="k", linewidth=0.5)
        append!(patches, [edgepatch])

    end

    localPatches = mplcollections.PatchCollection(patches, cmap=cmap, match_original=true)

    ax.add_collection(localPatches)

    ax.scatter(mpasOcean.xVertex, mpasOcean.yVertex, s=dotsize, c=phi, cmap=cmap)

    close()
    return fig, ax, cbar, localPatches
end
