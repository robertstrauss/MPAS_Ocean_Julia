using DelimitedFiles
import Base.Threads
using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
root = 0
commsize = MPI.Comm_size(comm)

nXChunks = sqrt(commsize)
if round(nXChunks) != nXChunks
	if rank == root
		println("ERROR: comm size $(commsize) not a square number (sqrt=$(nXChunks))")
	end
	throw()
end
nXChunks = Int(round(nXChunks))

if rank == root
	println("running with $(commsize) processes and $(Threads.nthreads()) threads per process")
end

include("mode_init/MPAS_Ocean.jl")
include("mode_init/MPAS_OceanHalos.jl")

if rank == root
	println("loading ocean from file...")
end

# fullOcean = MPAS_Ocean("MPAS_O_Shallow_Water/Mesh+Initial_Condition+Registry_Files/Periodic", "base_mesh.nc", "mesh.nc", periodicity="Periodic")

fullOcean = MPAS_Ocean("/global/cscratch1/sd/mpeterse/runs/210930_coastal_kelvin_wave_scaling2_1000x1000/ocean/verification/coastal_kelvin_wave/scaling1/initial_state/", "initial_state.nc", "initial_state.nc", periodicity="NonPeriodic_x")
# fullOcean = MPAS_Ocean("/global/homes/r/rstrauss/repos/MPAS_Ocean_Julia/MPAS_O_Shallow_Water/MPAS_O_Shallow_Water_Mesh_Generation/CoastalKelvinWaveMesh",
#                 "culled_mesh.nc", "mesh.nc", periodicity="NonPeriodic_x")

halowidth = 3

if rank == root
	println("loaded ocean from file with $(fullOcean.nCells) cells")
end

if rank == root
	println("dividing ocean into $(nXChunks)x$(nXChunks) among $(commsize) processes")
end
cellsInChunk, edgesInChunk, verticesInChunk, cellsFromChunk, cellsToChunk = divide_ocean(fullOcean, halowidth, nXChunks, nXChunks)#; iChunk = rank+1)

open("/global/cscratch1/sd/rstrauss/runs/fortrancomparison4096/divided1000x1000ocean_$(rank).txt", "w") do io
	writedlm(io, [cellsInChunk, edgesInChunk, verticesInChunk, cellsFromChunk, cellsToChunk])
end

myCells = cellsInChunk[rank+1] # cellsedgesvertices[1]
myEdges = edgesInChunk[rank+1] # cellsedgesvertices[2]
myVertices = verticesInChunk[rank+1] # cellsedgesvertices[3]

myOcean = mpas_subset(fullOcean, myCells, myEdges, myVertices)

MPI.Barrier(comm)

if rank == root
	println("ocean distributed between ranks.")
end

############### Ocean is now distributed between nodes. Let's do some simulation ####################

# set up initial condition

# include("mode_init/initial_conditions.jl")

mpasOcean = myOcean
meanCoriolisParameterf = sum(mpasOcean.fEdge) / length(mpasOcean.fEdge)
meanFluidThicknessH = sum(mpasOcean.bottomDepth)/length(mpasOcean.bottomDepth)
c = sqrt(mpasOcean.gravity*meanFluidThicknessH)
rossbyRadiusR = c/meanCoriolisParameterf

lYedge = maximum(fullOcean.yEdge) - minimum(fullOcean.yEdge)

function lateralProfile(y)
    return 1e-3*cos(y/lYedge * 4 * pi)
end

function kelvinWaveExactNormalVelocity(mpasOcean, iEdge, t=0)
    v = sqrt(mpasOcean.gravity*meanFluidThicknessH) * lateralProfile(mpasOcean.yEdge[iEdge] .+ c*t) * exp(-mpasOcean.xEdge[iEdge]/rossbyRadiusR)
    return v*sin(mpasOcean.angleEdge[iEdge])
end

function kelvinWaveExactSSH(mpasOcean, iCell, t=0)
    return - meanFluidThicknessH * lateralProfile(mpasOcean.yCell[iCell] .+ c*t) * exp(-mpasOcean.xCell[iCell]/rossbyRadiusR)
end

function kelvinWaveExactSolution!(mpasOcean, t=0)
    Threads.@threads for iCell in 1:mpasOcean.nCells
        mpasOcean.sshCurrent[iCell] = kelvinWaveExactSSH(mpasOcean, iCell, t)
    end

    Threads.@threads for iEdge in 1:mpasOcean.nEdges
        mpasOcean.normalVelocityCurrent[iEdge] = kelvinWaveExactNormalVelocity(mpasOcean, iEdge, t)
    end
end

function boundaryCondition!(mpasOcean, t)
    Threads.@threads for iEdge in 1:mpasOcean.nEdges
        if mpasOcean.boundaryEdge[iEdge] == 1.0
            mpasOcean.normalVelocityCurrent[iEdge] = kelvinWaveExactNormalVelocity(mpasOcean, iEdge, t)
        end
    end

end

# gaussianInit!(myOcean)

kelvinWaveExactSolution!(myOcean, 0)



MPI.Barrier(comm)

if rank == root
	println("initial condition set.")
end

# simulate for a while

include("mode_forward/time_steppers.jl")

nFrames = 100

# only use this core's subset, not entire ocean
mpasOcean = myOcean

# sshOverTime = zeros(nFrames+1, mpasOcean.nCells)
# sshOverTime[1,:] = mpasOcean.sshCurrent


if rank == root
	println("now simulating for $(nFrames) timesteps, updating halos every $(1) steps.")
end


@time begin
	exacttime = 0
	for f in 1:nFrames
		global exacttime
	# simulate until the halo areas are all invalid and need to be updated
	# for h in 1:6
	# forward_backward_step_threads!(mpasOcean)
        calculate_normal_velocity_tendency_threads!(mpasOcean)
    
        update_normal_velocity_by_tendency_threads!(mpasOcean)
    
	boundaryCondition!(mpasOcean, exacttime)

        calculate_ssh_tendency_threads!(mpasOcean)
    
        update_ssh_by_tendency_threads!(mpasOcean)
	# end
	
	exacttime += mpasOcean.dt

	### request cells in my halo from chunks with those cells
	halobufferssh = [] # temporarily stores new halo ssh
	halobuffernv = [] # temporarily stores new halo normal velocity
	recreqs = []
	for (srcchunk, localcells) in cellsFromChunk[rank+1]	
		newhalossh = Array{eltype(mpasOcean.sshCurrent)}(undef, length(localcells))
		append!(halobufferssh, [newhalossh])
		reqssh = MPI.Irecv!(newhalossh, srcchunk-1, 0, comm) # tag 0 for ssh
		append!(recreqs, [reqssh])

		localedges = collect(Set(mpasOcean.edgesOnCell[:,localcells]))
		newhalonv = Array{eltype(mpasOcean.normalVelocityCurrent)}(undef, length(localedges))
		append!(halobuffernv, [newhalonv])
		reqnv = MPI.Irecv!(newhalonv, srcchunk-1, 1, comm) # tag 1 for norm vel
		append!(recreqs, [reqnv])
	end
	
	MPI.Barrier(comm)
	### send cells in main non-halo area to chunks that need them for their halo
	sendreqs = []
	for (dstchunk, localcells) in cellsToChunk[rank+1]
		reqssh = MPI.Isend(mpasOcean.sshCurrent[localcells], dstchunk-1, 0, comm)
		append!(sendreqs, [reqssh])
		
		localedges = collect(Set(mpasOcean.edgesOnCell[:,localcells])) # Set to remove duplicates
		order = sortperm(myEdges[localedges])
		reqnv = MPI.Isend(mpasOcean.normalVelocityCurrent[localedges[order]], dstchunk-1, 1, comm)
		append!(sendreqs, [reqnv])
	end
	
	### copy the recieved data into the ocean's halo
	if rank == root
	# 	println("halo buffer before: ", halobuffernv[1][1:10])
	end
	MPI.Barrier(comm)
	MPI.Waitall!([recreqs..., sendreqs...])
	if rank == root
	# 	println("halo buffer after: ", halobuffernv[1][1:10])
	end
	MPI.Barrier(comm)
	for (i, (_, localcells)) in enumerate(cellsFromChunk[rank+1])
		mpasOcean.sshCurrent[localcells] = halobufferssh[i]
		localedges = collect(Set(mpasOcean.edgesOnCell[:,localcells]))
		order = sortperm(myEdges[localedges])
		mpasOcean.normalVelocityCurrent[localedges[order]] = halobuffernv[i]
	end

	# sshOverTime[f+1,:] = mpasOcean.sshCurrent

# 	if rank == root
# 		println("iteration $f of $nFrames complete")
# 	end
end end

if rank == root
	println("cells per proc: ", length(myOcean.sshCurrent))
end
# if rank == root
#	println("simulation complete, writing output arrays")
# end

# using DelimitedFiles

# open("out/dist_test_1_ssh_rank_$rank.txt", "w") do io
# 	writedlm(io, sshOverTime)
# end

# println("$rank: file written")
