include("mode_forward/time_steppers.jl")
include("mode_init/MPAS_Ocean.jl")
include("visualization.jl")

using DelimitedFiles

nCellsX = 100

mpasOcean = MPAS_Ocean(CODE_ROOT * "MPAS_O_Shallow_Water/MPAS_O_Shallow_Water_Mesh_Generation/CoastalKelvinWaveMesh/ConvergenceStudyMeshes",
                    "culled_mesh_$(nCellsX).nc",
                    "mesh_$(nCellsX).nc", periodicity="NonPeriodic_x", nvlevel=1)


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

nSamples = 10
nSteps = 24

sampletimes = zeros(Float64, nSamples)

exacttime = 0; kelvinWaveExactSolution!(mpasOcean, exacttime)

for sample in 1:nSamples
    sampletimes[sample] = @elapsed begin
        for i in 1:nSteps
            calculate_normal_velocity_tendency!(mpasOcean)

            boundaryCondition!(mpasOcean, exacttime)

            calculate_ssh_tendency!(mpasOcean)
            update_ssh_by_tendency!(mpasOcean)

            exacttime += mpasOcean.dt
        end
    end
end


fpath = CODE_ROOT * "output/serialCPU_timing/coastal_kelvinwave/steps_$nSteps/nCellsX_$nCellsX/"
mkpath(fpath)
fname = "$fpath$(Dates.now()).txt"
open(fname, "w") do io
    writedlm(io, sampletimes)
end
println("saved to $fname")

