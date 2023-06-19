import MPI


include("../mode_init/MPAS_Ocean.jl")


function update_halos!(comm, mpasOcean, haloMesh) # cellsFromMyRank, cellsToMyChunk, edgesFromMyChunk, edgesToMyChunk,)
	### request ssh & norm vel of cells in my halo from ranks with those cells
	recreqs = Vector{MPI.Request}()

	for (srcrank, localcells) in haloMesh.cellsToMyRankFrom
		reqthick = MPI.Irecv!(haloMesh.haloBuffersCells[srcrank], srcrank-1, 0, comm) # tag 0 for ~~ssh~~ thickness
		append!(recreqs, [reqthick])
	end
	for (srcrank, localedges) in haloMesh.edgesToMyRankFrom
		reqnv = MPI.Irecv!(haloMesh.haloBuffersEdges[srcrank], srcrank-1, 1, comm) # tag 1 for norm vel
		append!(recreqs, [reqnv])
	end


	MPI.Barrier(comm)

	### send cells in main non-halo area to ranks that need them for their halo
	sendreqs = Vector{MPI.Request}()

	for (dstrank, localcells) in haloMesh.cellsFromMyRankTo
		reqthick = MPI.Isend(mpasOcean.layerThickness[:,localcells], dstrank-1, 0, comm)
		append!(sendreqs, [reqthick])
	end
	for (dstrank, localedges) in haloMesh.edgesFromMyRankTo
		reqnv = MPI.Isend(mpasOcean.normalVelocityCurrent[:,localedges], dstrank-1, 1, comm)
		append!(sendreqs, [reqnv])
	end

# 	for (i, req) in enumerate(vcat(recreqs, sendreqs))
# 		# if MPI.Comm_rank(comm) == 0
# 		println("(rank $(MPI.Comm_rank(comm))) request $i status: $(MPI.Test!(req))")
# 		# end
# 	end

	MPI.Waitall!(vcat(recreqs, sendreqs))

	# copy data from buffers into ocean object
	for (srcrank, localcells) in haloMesh.cellsToMyRankFrom
		mpasOcean.layerThickness[:,localcells] = haloMesh.haloBuffersCells[srcrank]
	end
	for (srcrank, localedges) in haloMesh.edgesToMyRankFrom
		mpasOcean.normalVelocityCurrent[:,localedges] = haloMesh.haloBuffersEdges[srcrank] 
	end
end
