import CUDA
using SparseArrays

include("fixAngleEdge.jl")

"""
	GPU-located copy of the prognostic fields and arrays used in computation to advance prog. fields.
	Does not contain everything MPAS_Ocean does, only these fields.
"""
mutable struct MPAS_Ocean_CUDA

	# cell centered fields
	sshTendency::CUDA.CuArray{Float64,1}
	sshCurrent::CUDA.CuArray{Float64,1}
	bottomDepth::CUDA.CuArray{Float64,1}
	areaCell::CUDA.CuArray{Float64,1}


	nEdgesOnCell::CUDA.CuArray{Int64,1}

	edgesOnCell::CUDA.CuArray{Int64,2}
	cellsOnCell::CUDA.CuArray{Int64,2}

	edgeSignOnCell::CUDA.CuArray{Int64,2}

	# edge centered fields
	fEdge::CUDA.CuArray{Float64,1}
	dcEdge::CUDA.CuArray{Float64,1}
	dvEdge::CUDA.CuArray{Float64,1}
	yEdge::CUDA.CuArray{Float64,1}
	xEdge::CUDA.CuArray{Float64,1}
	angleEdge::CUDA.CuArray{Float64,1}


	normalVelocityTendency::CUDA.CuArray{Float64,2}
	normalVelocityCurrent::CUDA.CuArray{Float64,2}
	weightsOnEdge::CUDA.CuArray{Float64,2}


	nEdgesOnEdge::CUDA.CuArray{Int64,1}

	cellsOnEdge::CUDA.CuArray{Int64,2}
	edgesOnEdge::CUDA.CuArray{Int64,2}

	boundaryEdge::CUDA.CuArray{Int64,2} # only 0 or 1, could be less exact datatype. But will that slow down more by causing casting?

	# scalars
	nCells::Int64
	nEdges::Int64

	gravity::Float64
	dt::Float64


	function MPAS_Ocean_CUDA(mpasOcean::MPAS_Ocean)# where Float64<:AbstractFloat

		mpasOceanCuda = new()

		# mainArrayNames = [ # arrays used in computation arvancing prog. vars
	    #     :sshTendency,
	    #     :sshCurrent,
	    #     :nEdgesOnCell,
	    #     :edgesOnCell,
	    #     :cellsOnCell,
	    #     :bottomDepth,
	    #     :areaCell,
	    #     :edgeSignOnCell,
	    #     :normalVelocityTendency,
	    #     :normalVelocityCurrent,
	    #     :cellsOnEdge,
	    #     :nEdgesOnEdge,
	    #     :edgesOnEdge,
	    #     :weightsOnEdge,
	    #     :fEdge,
	    #     :dcEdge,
	    #     :dvEdge,
	    #     :boundaryEdge,
	    #     :yEdge,
	    #     :xEdge,
	    #     :angleEdge
	    # ]

		# scalarNames = [
		# 	:nCells,
		# 	:nEdges,
		# 	:gravity,
		# 	:dt
		# ]

	    for (field, typ) in zip(fieldnames(MPAS_Ocean_CUDA), MPAS_Ocean_CUDA.types)
	        setfield!(mpasOceanCuda, field, convert(typ, getfield(mpasOcean, field)))
	    end

		return mpasOceanCuda
	end
end



### GPU time-advancing simulation






function forward_backward_step_cuda!(mpasOcean::MPAS_Ocean_CUDA)
    calculate_normal_velocity_tendency_cuda!(mpasOcean)

    update_normal_velocity_by_tendency_cuda!(mpasOcean)

    calculate_ssh_tendency_cuda!(mpasOcean)

    update_ssh_by_tendency_cuda!(mpasOcean)
end



function forward_euler_step_cuda!(mpasOcean::MPAS_Ocean_CUDA)
    calculate_normal_velocity_tendency_cuda!(mpasOcean)

    calculate_ssh_tendency_cuda!(mpasOcean)

    update_normal_velocity_by_tendency_cuda!(mpasOcean)

    update_ssh_by_tendency_cuda!(mpasOcean)
end




function calculate_normal_velocity_tendency_cuda!(mpasOcean::MPAS_Ocean_CUDA)
    CUDA.@cuda blocks=cld(mpasOcean.nEdges, 1024) threads=1024 maxregs=64 calculate_normal_velocity_tendency_cuda_kernel!(
                                                                        mpasOcean.nEdges,
                                                                        mpasOcean.normalVelocityTendency,
                                                                        mpasOcean.normalVelocityCurrent,
                                                                        mpasOcean.sshCurrent,
                                                                        mpasOcean.cellsOnEdge,
                                                                        mpasOcean.nEdgesOnEdge,
                                                                        mpasOcean.edgesOnEdge,
                                                                        mpasOcean.weightsOnEdge,
                                                                        mpasOcean.fEdge,
                                                                        mpasOcean.dcEdge,
                                                                        mpasOcean.gravity)
end

function calculate_normal_velocity_tendency_cuda_kernel!(nEdges,
                                                          normalVelocityTendency,
                                                          normalVelocity,
                                                          ssh,
                                                          cellsOnEdge,
                                                          nEdgesOnEdge,
                                                          edgesOnEdge,
                                                          weightsOnEdge,
                                                          fEdge,
                                                          dcEdge,
                                                          gravity)

    iEdge = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
    if iEdge <= nEdges


        # gravity term: take gradient of ssh across edge
        cell1 = cellsOnEdge[1,iEdge]
        cell2 = cellsOnEdge[2,iEdge]

        if cell1 != 0 && cell2 != 0
            normalVelocityTendency[iEdge] = gravity * ( ssh[cellsOnEdge[1,iEdge]] - ssh[cellsOnEdge[2,iEdge]] ) / dcEdge[iEdge]
        end

        # coriolis term: combine norm. vel. of surrounding edges to approx. tangential vel.
        for i = 1:nEdgesOnEdge[iEdge]
            eoe = edgesOnEdge[i, iEdge]

            if eoe != 0
                normalVelocityTendency[iEdge] += weightsOnEdge[i,iEdge] * normalVelocity[eoe] * fEdge[eoe]
            end
        end


    end
    return
end







function update_normal_velocity_by_tendency_cuda!(mpasOcean::MPAS_Ocean_CUDA)
    CUDA.@cuda blocks=cld(mpasOcean.nEdges, 1024) threads=1024 maxregs=64 update_normal_velocity_by_tendency_cuda_kernel!(
                                                                        mpasOcean.nEdges,
                                                                        mpasOcean.normalVelocityCurrent,
                                                                        mpasOcean.dt,
                                                                        mpasOcean.normalVelocityTendency)
end

function update_normal_velocity_by_tendency_cuda_kernel!(nEdges,
                                                         normalVelocityCurrent,
                                                         dt,
                                                         normalVelocityTendency)
    iEdge = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
    if iEdge <= nEdges

        normalVelocityCurrent[iEdge] += dt * normalVelocityTendency[iEdge]

    end
    return
end








function calculate_ssh_tendency_cuda!(mpasOcean::MPAS_Ocean_CUDA)
    CUDA.@cuda blocks=cld(mpasOcean.nCells, 1024) threads=1024 maxregs=64 calculate_ssh_tendency_cuda_kernel!(
                                                                        mpasOcean.nCells,
                                                                        mpasOcean.sshTendency,
                                                                        mpasOcean.sshCurrent,
                                                                        mpasOcean.normalVelocityCurrent,
                                                                        mpasOcean.bottomDepth,
                                                                        mpasOcean.nEdgesOnCell,
                                                                        mpasOcean.edgesOnCell,
                                                                        mpasOcean.cellsOnCell,
                                                                        mpasOcean.areaCell,
                                                                        mpasOcean.edgeSignOnCell,
                                                                        mpasOcean.dvEdge)
end

function calculate_ssh_tendency_cuda_kernel!(nCells,
                                              sshTendency,
                                              ssh,
                                              normalVelocity,
                                              bottomDepth,
                                              nEdgesOnCell,
                                              edgesOnCell,
                                              cellsOnCell,
                                              areaCell,
                                              edgeSignOnCell,
                                              dvEdge)

    iCell = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
    if iCell <= nCells

        sshTendency[iCell] = 0

        # sum flux through each edge of cell
        for i in 1:nEdgesOnCell[iCell]
            edgeID = edgesOnCell[i,iCell]
            neighborCellID = cellsOnCell[i,iCell]

            if neighborCellID != 0
                depth = ( bottomDepth[neighborCellID] + bottomDepth[iCell] ) / 2.0
            else
                depth = bottomDepth[iCell]
            end

            flux = depth * normalVelocity[edgeID]
            sshTendency[iCell] += flux * edgeSignOnCell[iCell,i] * dvEdge[edgeID] / areaCell[iCell]
        end

    end
    return
end


function update_ssh_by_tendency_cuda!(mpasOcean::MPAS_Ocean_CUDA)
    CUDA.@cuda blocks=cld(mpasOcean.nCells, 1024) threads=1024 maxregs=64 update_ssh_by_tendency_cuda_kernel!(
                                                                        mpasOcean.nCells,
                                                                        mpasOcean.sshCurrent,
                                                                        mpasOcean.dt,
                                                                        mpasOcean.sshTendency)
end

function update_ssh_by_tendency_cuda_kernel!(nCells,
                                             sshCurrent,
                                             dt,
                                             sshTendency)
    iCell = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
    if iCell <= nCells

        sshCurrent[iCell] += dt * sshTendency[iCell]

    end
    return
end
