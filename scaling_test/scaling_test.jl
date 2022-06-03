using MPI
MPI.Init()

comm = MPI.COMM_WORLD
worldrank = MPI.Comm_rank(comm) + 1
root = 1
commsize = MPI.Comm_size(comm)




function rootprint(str)
	if worldrank == root
		println(str)
	end
end


using DelimitedFiles
using DataFrames
using CSV
using LinearAlgebra
using BenchmarkTools
import Plots
import Dates


CODE_ROOT = pwd() * "/"#../"
include(CODE_ROOT * "mode_init/MPAS_Ocean.jl")

include(CODE_ROOT * "mode_init/MPAS_OceanHalos.jl")

include(CODE_ROOT * "mode_init/initial_conditions.jl")
include(CODE_ROOT * "mode_init/exactsolutions.jl")

include(CODE_ROOT * "mode_forward/update_halos.jl")
include(CODE_ROOT * "mode_forward/time_steppers.jl")




function rectangularfactor(n)
	x = Int(round(sqrt(n)))
	if x^2 == n
		return x, x
	end
	while x > 0 && n/x !== round(n/x)
		x -= 1 end
	return x, Int(round(n/x))
end



function globalToLocal(globalcells, mycells)
	localcells = findall(iCell -> iCell in globalcells, mycells)
	return localcells
end



function runtests(proccounts; nsamples=6, nCellsX=64, halowidth=5, ncycles=2, nvlevels=100)
	global fname, partitiondir;

	rootprint("running tests for $proccounts counts of processors with $nsamples samples per test of $(halowidth*ncycles) steps each")


	fullOcean = MPAS_Ocean(CODE_ROOT * "MPAS_O_Shallow_Water/MPAS_O_Shallow_Water_Mesh_Generation/CoastalKelvinWaveMesh/ConvergenceStudyMeshes",
			       "culled_mesh_$nCellsX.nc", "mesh_$nCellsX.nc", periodicity="NonPeriodic_x", nvlevels=nvlevels)
	lYedge = maximum(fullOcean.yEdge) - minimum(fullOcean.yEdge)
	lateralProfilePeriodic(y) = 1e-3*cos(y/lYedge * 4 * pi)
	kelvinWaveExactNV, kelvinWaveExactSSH, kelvinWaveExactSolution!, boundaryCondition! = kelvinWaveGenerator(fullOcean, lateralProfilePeriodic)

	df = DataFrame([Int64[], collect([Float64[] for i in 1:nsamples])...,     collect([Float64[] for i in 1:nsamples])...,     Float64[], Float64[]],
		       ["procs", collect(["sim_time$i" for i in 1:nsamples])...,  collect(["mpi_time$i" for i in 1:nsamples])...,  "max_error", "l2_error"])

	for nprocs in proccounts
		splitcomm = MPI.Comm_split(comm, Int(worldrank>nprocs), worldrank) # only use the necessary ranks for this trial
		if worldrank <= nprocs


			partitionfile = partitiondir * "graph.info.part.$nprocs"
			rootprint("partitioning $nCellsX x $nCellsX x $(fullOcean.nVertLevels) ocean among $nprocs processes with $(partitionfile)")
			partitions = dropdims(readdlm(partitionfile), dims=2)
			mycells = findall(proc -> proc == worldrank-1, partitions)

			haloCells = grow_halo(fullOcean, mycells, 5)
			# rootprint("\t halo $haloCells")

			cellsFromChunk = Dict{Int64, Array}()
			for iCell in haloCells
				chunk = partitions[iCell] + 1
				if chunk in keys(cellsFromChunk)
					push!(cellsFromChunk[chunk], iCell)
				else
					cellsFromChunk[chunk] = [iCell]
				end
			end

			# tell neighbors which cells they should send me
			cellsToChunk = Dict{Int64, Array}()
			tag = 2
			sendreqs = []; recvreqs = []
			for (chunk, cells) in cellsFromChunk
				push!(sendreqs, MPI.Isend(cells, chunk-1, tag, splitcomm))

				numcells = MPI.Get_count(MPI.Probe(chunk-1, tag, splitcomm), Int64)
				cellsToChunk[chunk] = Array{Int64,1}(undef, numcells)
				push!(recvreqs, MPI.Irecv!(cellsToChunk[chunk], chunk-1, tag, splitcomm))
			end
			MPI.Waitall!([sendreqs..., recvreqs...])


			# rootprint("doing rect factor")
			# nx, ny = rectangularfactor(nprocs)
			# rootprint("dividing  $nCellsX x $nCellsX x $(fullOcean.nVertLevels) ocean among $nprocs processes with $nx x $ny grid")
			# cellsInChunk, edgesInChunk, verticesInChunk, cellsFromChunk, cellsToChunk = divide_ocean(fullOcean, halowidth, nx, ny)


			# myCells = cellsInChunk[worldrank] # cellsedgesvertices[1]
			# myEdges = edgesInChunk[worldrank] # cellsedgesvertices[2]
			# myVertices = verticesInChunk[worldrank] # cellsedgesvertices[3]

			rootprint("subsetting ocean")

			myCells     = union(mycells, haloCells)
	        myEdges     = collect(Set(fullOcean.edgesOnCell[:,myCells]))
	        myVertices  = collect(Set(fullOcean.verticesOnCell[:,myCells]))
			mpasOcean   = mpas_subset(fullOcean, myCells, myEdges, myVertices)

			cellsToMyChunk = Dict{Int64, Array}()
			for (chunk, cells) in cellsToChunk
				localcells = globalToLocal(cells, myCells)
				if length(localcells) > 0
					order = sortperm(myCells[localcells]) # order according to global index, so received in correct place
					cellsToMyChunk[chunk] = localcells[order]
				end
			end
			cellsFromMyChunk = Dict{Int64, Array}()
			for (chunk, cells) in cellsFromChunk
				localcells = globalToLocal(cells, myCells)
				if length(localcells) > 0
					order = sortperm(myCells[localcells]) # order according to global index, so received in correct place
					cellsFromMyChunk[chunk] = localcells[order]
				end
			end
			edgesToMyChunk = Dict{Int64, Array}()
			for (chunk, localcells) in cellsToMyChunk
				localedges = collect(Set(mpasOcean.edgesOnCell[:,localcells])) # Set() to remove duplicates
				if length(localedges) > 0
					order = sortperm(myEdges[localedges])  # order according to global index, so received in correct place
					edgesToMyChunk[chunk] = localedges[order]
				end
			end
			edgesFromMyChunk = Dict{Int64, Array}()
			for (chunk, localcells) in cellsFromMyChunk
				localedges = collect(Set(mpasOcean.edgesOnCell[:,localcells])) # Set() to remove duplicates
				if length(localedges) > 0
					order = sortperm(myEdges[localedges])  # order according to global index, so received in correct place
					edgesFromMyChunk[chunk] = localedges[order]
				end
			end

			# if worldrank == 1
				# println(" \t ($worldrank) from: $([chunk => length(cells) for (chunk, cells) in cellsFromMyChunk])")
				# println(" \t ($worldrank) to: $([chunk => length(cells) for (chunk, cells) in cellsToMyChunk])")
			# end

			rootprint("now simulating for $(halowidth*ncycles) steps, communicating every $halowidth steps")

			sampletimessim = zeros(nsamples)
			sampletimesmpi = zeros(nsamples)
			exacttime=0; kelvinWaveExactSolution!(mpasOcean, exacttime)
			for jsample in 1:nsamples
				for f in 1:ncycles
					# simulate
					sampletimessim[jsample] += @elapsed begin
						for h in 1:halowidth
							calculate_normal_velocity_tendency!(mpasOcean)
							update_normal_velocity_by_tendency!(mpasOcean)

							boundaryCondition!(mpasOcean, exacttime)

							calculate_ssh_tendency!(mpasOcean)
							update_ssh_by_tendency!(mpasOcean)
						end
						exacttime += halowidth*mpasOcean.dt
					end
					sampletimesmpi[jsample] += @elapsed begin
						update_halos!(splitcomm, mpasOcean, cellsFromMyChunk, cellsToMyChunk, edgesFromMyChunk, edgesToMyChunk)
					end
				end
			end
			#=
			bench = @benchmark begin
					# simulate
					for f in 1:$ncycles
						for h in 1:$halowidth
							$calculate_normal_velocity_tendency!($mpasOcean)
							$update_normal_velocity_by_tendency!($mpasOcean)

							$boundaryCondition!($mpasOcean, exacttime)

							$calculate_ssh_tendency!($mpasOcean)
							$update_ssh_by_tendency!($mpasOcean)
						end
						exacttime += $halowidth*$mpasOcean.dt

						update_halos!($splitcomm, $worldrank, $mpasOcean, $cellsFromChunk, $cellsToChunk, $myCells, $myEdges, $myVertices)
					end
				end samples=nsamples setup=(exacttime=0; $kelvinWaveExactSolution!($mpasOcean, exacttime))
			sampletimes = bench.times=#


			rootprint("setting up exact sol")
			exactssh = zero(mpasOcean.sshCurrent)
			for iCell in 1:mpasOcean.nCells
				exactssh[iCell,:] .= kelvinWaveExactSSH(mpasOcean, iCell, halowidth*ncycles*mpasOcean.dt)
			end

			rootprint("calculating error")
			difference = mpasOcean.sshCurrent .- exactssh
			MaxErrorNorm = norm(difference, Inf)
			L2ErrorNorm = norm(difference/sqrt(float(mpasOcean.nCells)))

			push!(df, [nprocs, sampletimessim..., sampletimesmpi..., MaxErrorNorm, L2ErrorNorm])

			rootprint(df[df.procs .<= nprocs, :])
			if worldrank == root
				CSV.write(fname, df)
			end
			rootprint("saved at $fname")

		end
		rootprint("freeing splitcomm")
		MPI.free(splitcomm)
		MPI.Barrier(comm)
	end



	return
end

halowidth = 1
ncycles = 10

proccounts = Int.(round.( 2 .^collect(1:log2(commsize)) ))
nCellsX = parse(Int64, ARGS[1])
nsamples = parse(Int64, ARGS[2])

partitiondir = CODE_ROOT * "scaling_test/graphparts/$(nCellsX)x$(nCellsX)/"
fpath = CODE_ROOT * "output/kelvinwave/resolution$(nCellsX)x$(nCellsX)/steps$(halowidth*ncycles)"
mkpath(fpath)
fname = fpath * "/$(Dates.now()).txt"
rootprint("results to: $fname")

runtests(proccounts; nCellsX=nCellsX, nsamples=nsamples, halowidth=halowidth, ncycles=ncycles)
