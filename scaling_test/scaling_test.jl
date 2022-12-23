using MPI
MPI.Init()

comm = MPI.COMM_WORLD
worldrank = MPI.Comm_rank(comm) + 1
root = 1
commsize = MPI.Comm_size(comm)


using DelimitedFiles
using DataFrames
using CSV
using LinearAlgebra
using BenchmarkTools
# import Plots
import Dates

using TimerOutputs


function rootprint(str)
	if worldrank == root
		println("$(Dates.Time(Dates.now())): ", str)
	end
end




CODE_ROOT = pwd() * "/"#../"
include(CODE_ROOT * "mode_init/MPAS_Ocean.jl")

include(CODE_ROOT * "mode_init/MPAS_OceanHalos.jl")

# include(CODE_ROOT * "mode_init/initial_conditions.jl")
include(CODE_ROOT * "mode_init/exactsolutions.jl")

include(CODE_ROOT * "mode_forward/update_halos.jl")
include(CODE_ROOT * "mode_forward/time_steppers.jl")



#=
function rectangularfactor(n)
	x = Int(round(sqrt(n)))
	if x^2 == n
		return x, x
	end
	while x > 0 && n/x !== round(n/x)
		x -= 1 end
	return x, Int(round(n/x))
end
=#


function globalToLocal(globalcells, mycells)
	localcells = findall(iCell -> iCell in globalcells, mycells)
	return localcells
end


function runtests(proccounts, fname, partitiondir; nsamples=6, nCellsX=64, halowidth=5, ncycles=2, nvlevels=100)

	meshpath = CODE_ROOT * "MPAS_O_Shallow_Water/MPAS_O_Shallow_Water_Mesh_Generation/CoastalKelvinWaveMesh/ConvergenceStudyMeshes"

	# kelvin wave initial condition and exact solution functions
	if worldrank == root
		fullOcean = MPAS_Ocean(meshpath, "culled_mesh_$nCellsX.nc", "mesh_$nCellsX.nc", periodicity="NonPeriodic_x", nvlevels=nvlevels)
		lYedge = maximum(fullOcean.yEdge) - minimum(fullOcean.yEdge)
		meanCoriolisParameterf = sum(fullOcean.fEdge) / length(fullOcean.fEdge)
		meanFluidThicknessH = sum(fullOcean.bottomDepth)/length(fullOcean.bottomDepth)
		c = sqrt(fullOcean.gravity*meanFluidThicknessH)
		rossbyRadiusR = c/meanCoriolisParameterf
		gravity = fullOcean.gravity
		dt = fullOcean.dt
	else
		lYedge = 0.0
		meanCoriolisParameterf = 0.0
		meanFluidThicknessH = 0.0
		c = 0.0
		rossbyRadiusR = 0.0
		gravity = 0.0
		dt = 0.0
	end
	MPI.Barrier(comm)
	lYedge, meanCoriolisParameterf, meanFluidThicknessH, c, rossbyRadiusR, gravity, dt = MPI.bcast([lYedge, meanCoriolisParameterf, meanFluidThicknessH, c, rossbyRadiusR, gravity, dt], root-1, comm)

	lateralProfilePeriodic(y) = 1e-6*cos(y/lYedge * 4 * pi)
	lateralProfile = lateralProfilePeriodic

	function kelvinWaveExactNormalVelocity(mpasOcean, iEdge, t=0)
		v = c * lateralProfile.(mpasOcean.yEdge[iEdge] .+ c*t) .* exp.(-mpasOcean.xEdge[iEdge]/rossbyRadiusR)
		return v .* sin.(mpasOcean.angleEdge[iEdge])
        end

	function kelvinWaveExactSSH(mpasOcean, iCell, t=0)
		return - meanFluidThicknessH * lateralProfile.(mpasOcean.yCell[iCell] .+ c*t) .* exp.(-mpasOcean.xCell[iCell]/rossbyRadiusR)
	end

	function kelvinWaveExactSolution!(mpasOcean, t=0)
		# just calculate exact solution once then copy it to lower layers
		sshperlayer = kelvinWaveExactSSH(mpasOcean, collect(1:mpasOcean.nCells), t) / mpasOcean.nVertLevels
		for k in 1:mpasOcean.nVertLevels
		    mpasOcean.layerThickness[k,:] .+= sshperlayer # add, don't replace, thickness contributes to level depth
		end

		nvperlayer = kelvinWaveExactNormalVelocity(mpasOcean, collect(1:mpasOcean.nEdges), t) / mpasOcean.nVertLevels
		for k in 1:mpasOcean.nVertLevels
		    mpasOcean.normalVelocityCurrent[k,:] .= nvperlayer
		end
	end

	function boundaryCondition!(mpasOcean, t)
		for iEdge in 1:mpasOcean.nEdges
		    if mpasOcean.boundaryEdge[iEdge] == 1.0
			mpasOcean.normalVelocityCurrent[:,iEdge] .= kelvinWaveExactNormalVelocity(mpasOcean, iEdge, t)/mpasOcean.nVertLevels
		    end
		end

	    end

#	timout = TimerOutput()

	rootprint("running tests for $proccounts counts of processors with $nsamples samples per test of $(halowidth*ncycles) steps each")
	rootprint("halos: $halowidth wide")

# 	@timeit timout "load ocean file" (fullOcean = MPAS_Ocean(CODE_ROOT * "MPAS_O_Shallow_Water/MPAS_O_Shallow_Water_Mesh_Generation/CoastalKelvinWaveMesh/ConvergenceStudyMeshes",
# 	"culled_mesh_$nCellsX.nc", "mesh_$nCellsX.nc", periodicity="NonPeriodic_x", nvlevels=nvlevels))
#	lYedge = maximum(fullOcean.yEdge) - minimum(fullOcean.yEdge)
#	lateralProfilePeriodic(y) = 1e-3*cos(y/lYedge * 4 * pi)
#	kelvinWaveExactNV, kelvinWaveExactSSH, kelvinWaveExactSolution!, boundaryCondition! = kelvinWaveGenerator(fullOcean, lateralProfilePeriodic)

	df = DataFrame([Int64[], collect([Float64[] for i in 1:nsamples])...,     collect([Float64[] for i in 1:nsamples])...,     Float64[], Float64[]],
		       ["procs", collect(["sim_time$i" for i in 1:nsamples])...,  collect(["mpi_time$i" for i in 1:nsamples])...,  "max_error", "l2_error"])

	cellsOnCell = Array{Int64,2}(replace(readdlm(partitiondir * "graph.info")[2:end,:], ""=>0))

	for nprocs in proccounts
		splitcomm = MPI.Comm_split(comm, Int(worldrank>nprocs), worldrank) # only use the necessary ranks for this trial
		if worldrank <= nprocs


			partitionfile = partitiondir * "graph.info.part.$nprocs"
			rootprint("partitioning $nCellsX x $nCellsX x $(nvlevels) ocean among $nprocs processes with $(partitionfile)")
			partitions = dropdims(readdlm(partitionfile), dims=2)
			mycells = findall(proc -> proc == worldrank-1, partitions)

#			@timeit timout "grow halo"
			(haloCells = grow_halo(cellsOnCell, mycells, halowidth))

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
			sendreqs = Array{MPI.Request,1}()
			recvreqs = Array{MPI.Request,1}()

			MPI.Barrier(splitcomm)
			for (chunk, cells) in cellsFromChunk
				push!(sendreqs, MPI.Isend(cells, chunk-1, tag, splitcomm))
			end
			MPI.Barrier(splitcomm)
			for (chunk, cells) in cellsFromChunk
				numcells = MPI.Get_count(MPI.Probe(chunk-1, tag, splitcomm), Int64)
				cellsToChunk[chunk] = Array{Int64,1}(undef, numcells)
				push!(recvreqs, MPI.Irecv!(cellsToChunk[chunk], chunk-1, tag, splitcomm))
			end
			MPI.Waitall!(vcat(sendreqs, recvreqs))

			rootprint("subsetting ocean")

#			@timeit timout "subset" begin
				myCells     = union(mycells, haloCells)
				rootprint("cell count $(length(myCells))")
				mpasOcean = MPAS_Ocean(meshpath, "culled_mesh_$nCellsX.nc", "mesh_$nCellsX.nc", periodicity="NonPeriodic_x", nvlevels=nvlevels, cells=myCells)
				mpasOcean.dt = dt #1e-2 * minimum(mpasOcean.dcEdge) / c
				rootprint("dt: $(mpasOcean.dt)")

				# rootprint("edgesoncell $(mpasOcean.edgesOnCell[1:10]")
#			end
				myEdges = collect(Set(mpasOcean.edgesOnCell[:,:]))

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

			rootprint("now simulating for $(halowidth*ncycles*nsamples) steps, communicating every $halowidth steps, timing $(nsamples) samples of $(halowidth*ncycles) steps")

			sampletimessim = zeros(nsamples)
			sampletimesmpi = zeros(nsamples)
			exacttime=0; kelvinWaveExactSolution!(mpasOcean, exacttime)
#			@timeit timout "simulation"
			for jsample in 1:nsamples
				for f in 1:ncycles
					# simulate
					sampletimessim[jsample] += @elapsed begin
						for h in 1:halowidth
							calculate_normal_velocity_tendency!(mpasOcean)
							update_normal_velocity_by_tendency!(mpasOcean)

							exacttime += mpasOcean.dt
							boundaryCondition!(mpasOcean, exacttime)

							calculate_thickness_tendency!(mpasOcean)
							update_thickness_by_tendency!(mpasOcean)
						end
					end
					sampletimesmpi[jsample] += @elapsed begin
						update_halos!(splitcomm, mpasOcean, cellsFromMyChunk, cellsToMyChunk, edgesFromMyChunk, edgesToMyChunk)
					end
				end
			end

			rootprint("generating exact sol")
			exactssh = zero(mpasOcean.sshCurrent)
			for iCell in 1:mpasOcean.nCells
				exactssh[iCell] = kelvinWaveExactSSH(mpasOcean, iCell, halowidth*ncycles*mpasOcean.dt)
			end

			rootprint("calculating error")
			mpasOcean.sshCurrent = dropdims(sum(mpasOcean.layerThickness, dims=1), dims=1) - mpasOcean.bottomDepth
			rootprint("num: $(mpasOcean.sshCurrent[1:10]), exact: $(exactssh[1:10])")
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
		MPI.free(splitcomm)
		MPI.Barrier(comm)

# 		if worldrank == root
#			show(timout)
		# end
	end



	return
end

halowidth = 1
ncycles = 10

nCellsX = parse(Int64, ARGS[1])
nsamples = parse(Int64, ARGS[2])
if length(ARGS) > 2
	proccounts = parse.(Int64, split(ARGS[3], ','))
else
	proccounts = Int.(round.( 2 .^collect(1:log2(commsize)) ))
end
nvlevels = 100

partitiondir = CODE_ROOT * "scaling_test/graphparts/$(nCellsX)x$(nCellsX)/"
fpath = CODE_ROOT * "output/kelvinwave/resolution$(nCellsX)x$(nCellsX)/procs$(commsize)/steps$(halowidth*ncycles)/nvlevels$(nvlevels)"
mkpath(fpath)
fname = fpath * "/$(Dates.now()).txt"
rootprint("results to: $fname")

runtests(proccounts, fname, partitiondir; nCellsX=nCellsX, nsamples=nsamples, halowidth=halowidth, ncycles=ncycles, nvlevels=nvlevels)
