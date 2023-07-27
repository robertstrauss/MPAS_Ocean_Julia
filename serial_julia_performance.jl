include("mode_forward/time_steppers.jl")
include("mode_forward/calculate_tendencies.jl")
include("mode_init/MPAS_Ocean.jl")
include("mode_init/exactsolutions.jl")
CODE_ROOT = pwd() * "/"
using DelimitedFiles
using Dates


function kelvinwaveserialperformance(fname, nCellsX, nSamples=10, nSteps=10)

    mpasOcean = MPAS_Ocean(CODE_ROOT * "MPAS_O_Shallow_Water/ConvergenceStudyMeshes/CoastalKelvinWave",
                        "culled_mesh_$(nCellsX)x$(nCellsX).nc",
                        "mesh_$(nCellsX)x$(nCellsX).nc", periodicity="NonPeriodic_x", nvlevels=1)


    meanCoriolisParameterf = sum(mpasOcean.fEdge) / length(mpasOcean.fEdge)
    meanFluidThicknessH = sum(mpasOcean.bottomDepth)/length(mpasOcean.bottomDepth)
    c = sqrt(mpasOcean.gravity*meanFluidThicknessH)
    rossbyRadiusR = c/meanCoriolisParameterf

    function lateralProfilePeriodic(y)
        return 1e-3*cos(y/mpasOcean.lY * 4 * pi)
    end

    period = mpasOcean.lY / (4*pi) /c

    lateralProfile = lateralProfilePeriodic

    kelvinWaveExactNormalVelocity, kelvinWaveExactSSH, kelvinWaveExactSolution!, boundaryCondition! = kelvinWaveGenerator(mpasOcean, lateralProfile)

    sampletimes = zeros(Float64, nSamples)
    
    println("now simulating")
    exacttime = 0; kelvinWaveExactSolution!(mpasOcean, exacttime)

    for sample in 1:nSamples
        sampletimes[sample] = @elapsed begin
            for i in 1:nSteps
                calculate_normal_velocity_tendency!(mpasOcean)

                boundaryCondition!(mpasOcean, exacttime)

                calculate_thickness_tendency!(mpasOcean)
                update_thickness_by_tendency!(mpasOcean)

                exacttime += mpasOcean.dt
            end
        end
    end

    println("times: $sampletimes")
    
    open(fname, "w") do io
        writedlm(io, sampletimes)
    end
    
    println("saved to $fname")
    
    return mpasOcean, exacttime
end


nCellsX = parse(Int64, ARGS[1])
nSamples = parse(Int64, ARGS[2])
nSteps = 10

fpath = CODE_ROOT * "output/serialCPU_timing/coastal_kelvinwave/steps_$nSteps/resolution_$(nCellsX)x$(nCellsX)/"
mkpath(fpath)
fname = "$fpath$(Dates.now()).txt"
println("output file: $fname")

kelvinwaveserialperformance(fname, nCellsX, nSamples, nSteps)
