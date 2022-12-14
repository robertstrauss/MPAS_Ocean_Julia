function kelvinWaveGenerator(mpasOcean, lateralProfile)
    meanCoriolisParameterf = sum(mpasOcean.fEdge) / length(mpasOcean.fEdge)
    meanFluidThicknessH = sum(mpasOcean.bottomDepth)/length(mpasOcean.bottomDepth)
    c = sqrt(mpasOcean.gravity*meanFluidThicknessH)
    rossbyRadiusR = c/meanCoriolisParameterf

    function kelvinWaveExactNormalVelocity(mpasOcean, iEdge, t=0)
        v = c * lateralProfile.(mpasOcean.yEdge[iEdge] .+ c*t) .* exp.(-mpasOcean.xEdge[iEdge]/rossbyRadiusR)
        return v .* sin.(mpasOcean.angleEdge[iEdge])
    end

    function kelvinWaveExactSSH(mpasOcean, iCell, t=0)
        return - meanFluidThicknessH * lateralProfile.(mpasOcean.yCell[iCell] .+ c*t) .* exp.(-mpasOcean.xCell[iCell]/rossbyRadiusR)
    end

    function kelvinWaveExactSolution!(mpasOcean, t=0)
        mpasOcean.sshCurrent[:] = kelvinWaveExactSSH(mpasOcean, collect(1:mpasOcean.nCells), t)

        mpasOcean.normalVelocityCurrent[:,1] = kelvinWaveExactNormalVelocity(mpasOcean, collect(1:mpasOcean.nEdges), t)
        for k in 2:mpasOcean.nVertLevels
            mpasOcean.normalVelocityCurrent[:,k] = mpasOcean.normalVelocityCurrent[:,1]
        end
    end

    function boundaryCondition!(mpasOcean, t)
        for iEdge in 1:mpasOcean.nEdges
            if mpasOcean.boundaryEdge[iEdge] == 1.0
                mpasOcean.normalVelocityCurrent[iEdge,:] .= kelvinWaveExactNormalVelocity(mpasOcean, iEdge, t)/mpasOcean.nVertLevels
            end
        end

    end

    return kelvinWaveExactNormalVelocity, kelvinWaveExactSSH, kelvinWaveExactSolution!, boundaryCondition!
end



function inertiaGravityWaveParams(mpasOcean, nX, nY, etaHat)
    f0 = sum(mpasOcean.fEdge) / length(mpasOcean.fEdge)
    
    meanFluidThickness = sum(mpasOcean.bottomDepth)/length(mpasOcean.bottomDepth)
    
    kX = nX * 2*pi / mpasOcean.lX

    kY = nY * 2*pi / mpasOcean.lY

    omega = sqrt(f0^2 + mpasOcean.gravity*meanFluidThickness*(kX^2 + kY^2))

    return etaHat, f0, kX, kY, omega
end
function inertiaGravityExactSolution!(mpasOcean::MPAS_Ocean, etaHat::Float64, f0::Float64, kX::Float64, kY::Float64, omega::Float64, t::Float64)
    for iCell in 1:mpasOcean.nCells
        x = mpasOcean.xCell[iCell]
        y = mpasOcean.yCell[iCell]
        
        mpasOcean.sshCurrent[iCell] = DetermineInertiaGravityWaveExactSurfaceElevation(etaHat,kX,kY,omega,x,y,t)
    end
    
    for iEdge in 1:mpasOcean.nEdges
        x = mpasOcean.xEdge[iEdge]
        y = mpasOcean.yEdge[iEdge]
        
        u = DetermineInertiaGravityWaveExactZonalVelocity(etaHat, f0, mpasOcean.gravity, kX, kY, omega, x, y, t)
        
        v = DetermineInertiaGravityWaveExactMeridionalVelocity(etaHat, f0, mpasOcean.gravity, kX, kY, omega, x, y, t)
        
        theta = mpasOcean.angleEdge[iEdge]
        
        mpasOcean.normalVelocityCurrent[iEdge] = u*cos(theta) + v*sin(theta)
    end
end
function inertiaGravityExactNormalVelocity(theta, etaHat,f0,g,kX,kY,omega,x,y,t)

    u = DetermineInertiaGravityWaveExactZonalVelocity(etaHat, f0, g, kX, kY, omega, x, y, t)

    v = DetermineInertiaGravityWaveExactMeridionalVelocity(etaHat, f0, g, kX, kY, omega, x, y, t)

    return u*cos(theta) + v*sin(theta)
end
function DetermineInertiaGravityWaveExactMeridionalVelocity(etaHat,f0,g,kX,kY,omega,x,y,time)
    v = etaHat*(g/(omega^2.0 - f0^2.0)*(omega*kY*cos(kX*x + kY*y - omega*time)
                                          + f0*kX*sin(kX*x + kY*y - omega*time)))
    return v
end
function DetermineInertiaGravityWaveExactZonalVelocity(etaHat,f0,g,kX,kY,omega,x,y,time)
    u = etaHat*(g/(omega^2.0 - f0^2.0)*(omega*kX*cos(kX*x + kY*y - omega*time)
                                          - f0*kY*sin(kX*x + kY*y - omega*time)))
    return u
end
function DetermineInertiaGravityWaveExactSurfaceElevation(etaHat,kX,kY,omega,x,y,time)
    eta = etaHat*cos(kX*x + kY*y - omega*time)
    return eta
end
