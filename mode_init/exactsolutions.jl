function kelvinWaveGenerator(mpasOcean, lateralProfile)
    meanCoriolisParameterf = sum(mpasOcean.fEdge) / length(mpasOcean.fEdge)
    meanFluidThicknessH = sum(mpasOcean.bottomDepth)/length(mpasOcean.bottomDepth)
    c = sqrt(mpasOcean.gravity*meanFluidThicknessH)
    rossbyRadiusR = c/meanCoriolisParameterf

    function kelvinWaveExactNormalVelocity(mpasOcean, iEdge, t=0)
        v = sqrt(mpasOcean.gravity*meanFluidThicknessH) * lateralProfile(mpasOcean.yEdge[iEdge] .+ c*t) * exp(-mpasOcean.xEdge[iEdge]/rossbyRadiusR)
        return v*sin(mpasOcean.angleEdge[iEdge])
    end

    function kelvinWaveExactSSH(mpasOcean, iCell, t=0)
        return - meanFluidThicknessH * lateralProfile(mpasOcean.yCell[iCell] .+ c*t) * exp(-mpasOcean.xCell[iCell]/rossbyRadiusR)
    end

    function kelvinWaveExactSolution!(mpasOcean, t=0)
        for iCell in 1:mpasOcean.nCells
            mpasOcean.sshCurrent[iCell] = kelvinWaveExactSSH(mpasOcean, iCell, t)
        end

        for iEdge in 1:mpasOcean.nEdges
            mpasOcean.normalVelocityCurrent[iEdge] = kelvinWaveExactNormalVelocity(mpasOcean, iEdge, t)
        end
    end

    function boundaryCondition!(mpasOcean, t)
        for iEdge in 1:mpasOcean.nEdges
            if mpasOcean.boundaryEdge[iEdge] == 1.0
                mpasOcean.normalVelocityCurrent[iEdge] = kelvinWaveExactNormalVelocity(mpasOcean, iEdge, t)
            end
        end

    end

    return kelvinWaveExactNormalVelocity, kelvinWaveExactSSH, kelvinWaveExactSolution!, boundaryCondition!
end
