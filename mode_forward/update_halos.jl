import MPI


include("../mode_init/MPAS_Ocean.jl")


function update_halos!(comm, mpasOcean, cellsFromMyChunk, cellsToMyChunk, myCells, myEdges)
	### request cells in my halo from chunks with those cells
	halobufferssh = [] # temporarily stores new halo ssh
	halobuffernv = [] # temporarily stores new halo normal velocity
	recreqs = Vector{MPI.Request}()
	for (srcchunk, localcells) in cellsFromMyChunk
		println(" \t \t receive req $srcchunk $(length(localcells))")
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
	sendreqs = Vector{MPI.Request}()
	for (dstchunk, localcells) in cellsToMyChunk
		println(" \t \t send    req $dstchunk $(length(localcells))")
		reqssh = MPI.Isend(mpasOcean.sshCurrent[localcells], dstchunk-1, 0, comm)
		append!(sendreqs, [reqssh])

		localedges = collect(Set(mpasOcean.edgesOnCell[:,localcells])) # Set to remove duplicates
		order = sortperm(myEdges[localedges])
		reqnv = MPI.Isend(mpasOcean.normalVelocityCurrent[localedges[order]], dstchunk-1, 1, comm)
		append!(sendreqs, [reqnv])
	end

	### copy the recieved data into the ocean's halo
	MPI.Barrier(comm)
	println(" rec $(length(recreqs))")
	println(" T send $(typeof(sendreqs)), $(eltype(sendreqs))")
	for (i, recreq) in enumerate(recreqs)
		println(i)
		flag, status = MPI.Test!(recreq)
		println("\t f s $flag \t $status")
	end
	MPI.Waitall!([recreqs..., sendreqs...])
	# MPI.Barrier(comm)
	for (i, (_, localcells)) in enumerate(cellsFromMyChunk)
		mpasOcean.sshCurrent[localcells] = halobufferssh[i]
		localedges = collect(Set(mpasOcean.edgesOnCell[:,localcells]))
		order = sortperm(myEdges[localedges])
		mpasOcean.normalVelocityCurrent[localedges[order]] = halobuffernv[i]
	end
end
