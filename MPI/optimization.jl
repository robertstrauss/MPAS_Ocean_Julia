
CODE_ROOT = pwd() * "/"#../"

include(CODE_ROOT * "mode_init/MPAS_Ocean.jl")
include(CODE_ROOT * "mode_init/MPAS_OceanHalos.jl")

include(CODE_ROOT * "mode_init/initial_conditions.jl")
include(CODE_ROOT * "mode_forward/time_steppers.jl")
include(CODE_ROOT * "mode_init/exactsolutions.jl")



F = Float64
I = Int64


# fullOcean = MPAS_Ocean("MPAS_O_Shallow_Water/Mesh+Initial_Condition+Registry_Files/Periodic", "base_mesh.nc", "mesh.nc", periodicity="Periodic")
nCellsX = 64::I
println("loading ocean from file")
fullOcean = MPAS_Ocean(CODE_ROOT * "MPAS_O_Shallow_Water/MPAS_O_Shallow_Water_Mesh_Generation/CoastalKelvinWaveMesh/ConvergenceStudyMeshes",
		       "culled_mesh_$nCellsX.nc", "mesh_$nCellsX.nc", periodicity="NonPeriodic_x")
halowidth = 5::I

println("dividing ocean")
cellsInChunk, edgesInChunk, verticesInChunk, cellsFromChunk, cellsToChunk = divide_ocean(fullOcean, halowidth, 2, 2)#; iChunk = rank)

rank = 1

myCells = cellsInChunk[rank] # cellsedgesvertices[1]
myEdges = edgesInChunk[rank] # cellsedgesvertices[2]
myVertices = verticesInChunk[rank] # cellsedgesvertices[3]

mpasOcean = mpas_subset(fullOcean, myCells, myEdges, myVertices)



# mpasOcean = fullOcean

lYedge = maximum(fullOcean.yEdge) - minimum(fullOcean.yEdge)
lateralProfilePeriodic(y) = 1e-3*cos(y/lYedge * 4 * pi)
kelvinWaveExactNV, kelvinWaveExactSSH, kelvinWaveExactSolution!, boundaryCondition! = kelvinWaveGenerator(mpasOcean, lateralProfilePeriodic)





using BenchmarkTools

exacttime = 0.0::F

@btime calculate_normal_velocity_tendency!(mpasOcean)
@btime update_normal_velocity_by_tendency!(mpasOcean)

@btime update_normal_velocity_by_tendency_loop!(mpasOcean)


@btime boundaryCondition!(mpasOcean, exacttime)

@btime calculate_ssh_tendency!(mpasOcean)
@btime update_ssh_by_tendency!(mpasOcean)

#
#


@code_warntype calculate_normal_velocity_tendency!(mpasOcean)
@code_warntype update_normal_velocity_by_tendency!(mpasOcean)

@code_warntype boundaryCondition!(mpasOcean, exacttime)

@code_warntype calculate_ssh_tendency!(mpasOcean)
