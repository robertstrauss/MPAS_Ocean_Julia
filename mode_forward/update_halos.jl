import MPI


include("../mode_init/MPAS_Ocean.jl")


function update_halos!(comm, mpasOcean, cellsFromMyChunk, cellsToMyChunk, edgesFromMyChunk, edgesToMyChunk,)
	### request ssh & norm vel of cells in my halo from chunks with those cells
	recreqs = Vector{MPI.Request}()

	halobufferssh = [] # temporarily stores new halo ssh
	for (srcchunk, localcells) in cellsFromMyChunk
		# println("ssh recv $(srcchunk-1) ($(length(localcells)))")
		newhalossh = Array{eltype(mpasOcean.sshCurrent)}(undef, length(localcells))
		append!(halobufferssh, [newhalossh])
		reqssh = MPI.Irecv!(newhalossh, srcchunk-1, 0, comm) # tag 0 for ssh
		append!(recreqs, [reqssh])
	end
	halobuffernv = [] # temporarily stores new halo normal velocity
	for (srcchunk, localedges) in edgesFromMyChunk
		newhalonv = Array{eltype(mpasOcean.normalVelocityCurrent)}(undef, length(localedges))
		append!(halobuffernv, [newhalonv])
		reqnv = MPI.Irecv!(newhalonv, srcchunk-1, 1, comm) # tag 1 for norm vel
		append!(recreqs, [reqnv])
	end


	MPI.Barrier(comm)

	### send cells in main non-halo area to chunks that need them for their halo
	sendreqs = Vector{MPI.Request}()

	for (dstchunk, localcells) in cellsToMyChunk
		reqssh = MPI.Isend(mpasOcean.sshCurrent[localcells], dstchunk-1, 0, comm)
		append!(sendreqs, [reqssh])
	end
	for (dstchunk, localedges) in edgesToMyChunk
		reqnv = MPI.Isend(mpasOcean.normalVelocityCurrent[localedges], dstchunk-1, 1, comm)
		append!(sendreqs, [reqnv])
	end

	### copy the recieved data into the ocean's halo
	MPI.Waitall!(vcat(recreqs, sendreqs))

	# copy data from buffers into ocean object
	for (i, (srcchunk, localcells)) in enumerate(cellsFromMyChunk)
		mpasOcean.sshCurrent[localcells] = halobufferssh[i]
	end
	for (i, (srcchunk, localedges)) in enumerate(edgesFromMyChunk)
		mpasOcean.normalVelocityCurrent[localedges] = halobuffernv[i]
	end
end
