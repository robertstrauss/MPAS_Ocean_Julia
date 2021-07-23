

function gaussianInit!(mpasOcean; sx=1/3, sy=1/3, mx=1/2, my=1/2)
    sigmax = mpasOcean.lX*sx
    sigmay = mpasOcean.lY*sy
    
    mux = mpasOcean.lX*mx
    muy = mpasOcean.lY*my
    
    for iCell in 1:mpasOcean.nCells
        mpasOcean.sshCurrent[iCell] = exp( - (mpasOcean.xCell[iCell] - mux)^2 / sigmax^2 - (mpasOcean.yCell[iCell] - muy)^2 / sigmay^2 )
    end
    
    mpasOcean.normalVelocityCurrent .= 0
    
    return nothing
end

function planeWaveInit!(mpasOcean; rkx = 1, rky = 1)
    kx = rkx * 2 * pi / mpasOcean.lX
    ky = rky * 2 * pi / mpasOcean.lY
    
    for iCell in 1:mpasOcean.nCells
        mpasOcean.sshCurrent[iCell] = cos( kx * mpasOcean.xCell[iCell] + ky * mpasOcean.yCell[iCell] )
    end
    
    mpasOcean.normalVelocityCurrent .= 0
    
    mpasOcean.sshTendency[:] .= mpasOcean.sshCurrent[:]
    mpasOcean.normalVelocityTendency[:] .= mpasOcean.normalVelocityCurrent[:]
    
    return nothing
end