import NCDatasets
import CUDA
using SparseArrays
# import KernelAbstractions

# include("Namelist.jl")
include("fixAngleEdge.jl")


# const I = Int64
# const F = Float64

function dimsortperm(A::AbstractMatrix; dims::Integer, rev::Bool = false)
    P = mapslices(x -> sortperm(x; rev = rev), A, dims = dims)
    if dims == 1
        for j = 1:size(P, 2)
            offset = (j - 1) * size(P, 1)
            for i = 1:size(P, 1)
                P[i, j] += offset
            end
        end
    else # if dims == 2
        for j = 1:size(P, 2)
            for i = 1:size(P, 1)
                P[i, j] = (P[i, j] - 1) * size(P, 1) + i
            end
        end
    end
    return P
end

mutable struct MPAS_Ocean#{Float64<:AbstractFloat}

#     boundaryConditions::Array{BoundaryCondition}

    # prognostic variables
    sshCurrent::Array{Float64,1}   #where Float64 <: AbstractFloat # sea surface height (cell-centered)
    sshTendency::Array{Float64,1}   #where Float64 <: AbstractFloat # tendency (cell-centered)

    normalVelocityCurrent::Array{Float64,2}   #where Float64 <: AbstractFloat # group velocity normal to mesh edges (edge-centered)
    normalVelocityTendency::Array{Float64,2}   #where Float64 <: AbstractFloat # tendency (edge-centered)


	bottomDepth::Array{Float64,1}   #where Float64 <: AbstractFloat # bathymetry (cell-centered)
    bottomDepthEdge::Array{Float64,1}   #where Float64 <: AbstractFloat # bathymetry (cell-centered)
    gravity::Float64

    dt::Float64 # timestep



    ### mesh information


    ## cell-centered arrays
    nCells::Int64  # number of cells
    cellsOnCell::Array{Int64,2}  # indices of neighbor cells (nCells, nEdgesOnCell[Int64])
    edgesOnCell::Array{Int64,2}  # indices of edges of cell (nCells, nEdgesOnCell[Int64])
    verticesOnCell::Array{Int64,2}  # indices of vertices of cell (nCells, nEdgesOnCell[Int64])
    kiteIndexOnCell::Array{Int64,2}  # index of kite shape on cell connecting vertex, midpoint of edge, center of cell, and midpoint of other edge (nCells, nEdgesOnCell[Int64])
    nEdgesOnCell::Array{Int64,1}  # number of edges a cell has (cell-centered)
    edgeSignOnCell::Array{Int8,2} # orientation of edge relative to cell (nCells, nEdgesOnCell[Int64])
    latCell::Array{Float64,1}   #where Float64 <: AbstractFloat # latitude of cell (cell-centered)
    lonCell::Array{Float64,1}   #where Float64 <: AbstractFloat # longitude of cell (cell-centered)
    xCell::Array{Float64,1}   #where Float64 <: AbstractFloat # x coordinate of cell (cell-centered)
    yCell::Array{Float64,1}   #where Float64 <: AbstractFloat # y coordinate of cell (cell-centered)
    areaCell::Array{Float64,1}   #where Float64 <: AbstractFloat # area of cell (cell-centered)
    fCell::Array{Float64,1}   #where Float64 <: AbstractFloat # coriolis parameter (cell-centered)
    maxLevelCell::Array{Int64,1}
    gridSpacing::Array{Float64,1}   #where Float64 <: AbstractFloat
    boundaryCell::Array{Int64,2}  # 0 for inner cells, 1 for boundary cells (cell-centered)
	cellMask::Array{Int64,2}


    ## edge-centered arrays
    nEdges::Int64
    cellsOnEdge::Array{Int64,2}  # indices of the two cells on this edge (nEdges, 2)
    edgesOnEdge::Array{Int64,2}
    verticesOnEdge::Array{Int64,2}
    nEdgesOnEdge::Array{Int64,1}
    xEdge::Array{Float64,1}   #where Float64 <: AbstractFloat
    yEdge::Array{Float64,1}   #where Float64 <: AbstractFloat
    dvEdge::Array{Float64,1}   #where Float64 <: AbstractFloat
    dcEdge::Array{Float64,1}   #where Float64 <: AbstractFloat
    fEdge::Array{Float64,1}   #where Float64 <: AbstractFloat # coriolis parameter
    angleEdge::Array{Float64,1}   #where Float64 <: AbstractFloat
    weightsOnEdge::Array{Float64,2}   #where Float64 <: AbstractFloat # coeffecients of norm vels of surrounding edges in linear combination to compute tangential velocity
    maxLevelEdgeTop::Array{Int64,1}
    maxLevelEdgeBot::Array{Int64,1}
    boundaryEdge::Array{Int64,2}
    edgeMask::Array{Int64,2}


    ## vertex-centered arrays
    nVertices::Int64
    latVertex::Array{Float64,1}   #where Float64 <: AbstractFloat
    lonVertex::Array{Float64,1}   #where Float64 <: AbstractFloat
    xVertex::Array{Float64,1}   #where Float64 <: AbstractFloat
    yVertex::Array{Float64,1}   #where Float64 <: AbstractFloat
    vertexDegree::Int64
    cellsOnVertex::Array{Int64}
    edgesOnVertex::Array{Int64}
    edgeSignOnVertex::Array{Int8}
    fVertex::Array{Float64}   #where Float64 <: AbstractFloat
    areaTriangle::Array{Float64}   #where Float64 <: AbstractFloat
    kiteAreasOnVertex::Array{Float64}   #where Float64 <: AbstractFloat
    maxLevelVertexTop::Array{Int64}
    maxLevelVertexBot::Array{Int64}
    boundaryVertex::Array{Int64,2}
    vertexMask::Array{Int64,2}



    nVertLevels::Int64
    layerThickness::Array{Float64,2}
    layerThicknessEdge::Array{Float64,2}


#     ExactSolutionParameters::Array{Float64}   #where Float64 <: AbstractFloat
#     myNamelist::Namelist

    gridSpacingMagnitude::Float64

    # maximum values of x and y (size of mesh)
    lX::Float64
    lY::Float64

    nNonPeriodicBoundaryCells::Int64
    nNonPeriodicBoundaryEdges::Int64
    nNonPeriodicBoundaryVertices::Int64

	function MPAS_Ocean(mesh_directory = "MPAS_O_Shallow_Water/Mesh+Initial_Condition+Registry_Files/Periodic",
						base_mesh_file_name = "base_mesh.nc",
						mesh_file_name = "mesh.nc";
						nvlevels = 1,
						periodicity = "Periodic",
						cells = Colon())# where Float64<:AbstractFloat

		# tmr = TimerOutput()

	# @timeit tmr "file" begin
		# @timeit tmr "ocean" begin
		mpasOcean = new()
		# end




		base_mesh_file = NCDatasets.Dataset("$mesh_directory/$base_mesh_file_name", "r", format=:netcdf4)
		mesh_file = NCDatasets.Dataset("$mesh_directory/$mesh_file_name", "r", format=:netcdf4)

		# Choose my_mesh_file_name to be either base_mesh_file_name or mesh_file_name
		my_mesh_file_name = mesh_file_name
		# Choose my_mesh_file to be either base_mesh_file or mesh_file
		my_mesh_file = mesh_file


		maxEdges = mesh_file.dim["maxEdges"]

		## defining the mesh
		nGlobalCells = mesh_file.dim["nCells"]
		nGlobalEdges = mesh_file.dim["nEdges"]
		nGlobalVertices = mesh_file.dim["nVertices"]
		mpasOcean.vertexDegree = mesh_file.dim["vertexDegree"]
		mpasOcean.nVertLevels = nvlevels

		# will be incremented in ocn_init_routines_compute_max_level
		mpasOcean.nNonPeriodicBoundaryEdges = 0
		mpasOcean.nNonPeriodicBoundaryVertices = 0
		mpasOcean.nNonPeriodicBoundaryCells = 0


		if cells == Colon()
			cells = 1:nGlobalCells
		end

		mpasOcean.edgesOnCell = my_mesh_file["edgesOnCell"][:]
		mpasOcean.verticesOnCell = my_mesh_file["verticesOnCell"][:]

		edges = collect(Set(mpasOcean.edgesOnCell[:,cells]))
		vertices = collect(Set(mpasOcean.verticesOnCell[:,cells]))
	# end
	# @timeit tmr "tolocal" begin


	    globtolocalcell = Dict{Int64,Int64}()
	    for (iLocal, iGlobal) in enumerate(cells)
		    globtolocalcell[iGlobal] = iLocal
	    end

	    globtolocaledge = Dict{Int64,Int64}()
	    for (iLocal, iGlobal) in enumerate(edges)
			# println("enumerate $iLocal $iGlobal")
		    globtolocaledge[iGlobal] = iLocal
	    end

	    globtolocalvertex = Dict{Int64,Int64}()
	    for (iLocal, iGlobal) in enumerate(vertices)
		    globtolocalvertex[iGlobal] = iLocal
	    end


		function globalToLocal!(globals, globmap)
			for i in 1:length(globals)
				if globals[i] in keys(globmap)
					globals[i] = globmap[globals[i]]
				else
					globals[i] = 0
				end
			end
			# sort!(globals, dims=1, rev=true)
			ordermap = dimsortperm(globals, dims=1, rev=true) # move "0" values to the end of each row, like a ragged array, like the code expects
			globals[:] = globals[ordermap]
			# # println("ordermap $ordermap")

			return ordermap # its now necesary to remember the way we reordered the indexes, and reorder related arrays the same way
		end

		mpasOcean.nCells = nCells = length(cells)
		mpasOcean.nVertices = nVertices = length(vertices)
		mpasOcean.nEdges = nEdges = length(edges)
# end
# @timeit tmr "mesh on cell" begin


		mpasOcean.cellsOnCell = my_mesh_file["cellsOnCell"][:][:,cells]#, globtolocalcell)
		mpasOcean.edgesOnCell = my_mesh_file["edgesOnCell"][:][:,cells]#, globtolocaledge)
		mpasOcean.verticesOnCell = my_mesh_file["verticesOnCell"][:][:,cells]#, globtolocalvertex)
# end
# @timeit tmr "mesh on edge" begin
		mpasOcean.cellsOnEdge = my_mesh_file["cellsOnEdge"][:][:,edges]# :: Array{Int64,2}#, globtolocalcell)
		# mpasOcean.nEdgesOnEdge = base_mesh_file["nEdgesOnEdge"][:][edges]
		mpasOcean.edgesOnEdge = my_mesh_file["edgesOnEdge"][:][:,edges]# :: Array{Int64,2}
		mpasOcean.verticesOnEdge = my_mesh_file["verticesOnEdge"][:][:,edges]# :: Array{Int64,2}
# end
# @timeit tmr "mesh on vertex" begin
		mpasOcean.cellsOnVertex = my_mesh_file["cellsOnVertex"][:][:,vertices]
		mpasOcean.edgesOnVertex = my_mesh_file["edgesOnVertex"][:][:,vertices]
		#

# end
# @timeit tmr "mut" begin
		# println("glob to local edge $globtolocaledge")
		# println("edges on cell pre $(mpasOcean.edgesOnCell[1:20])")

		globalToLocal!(mpasOcean.edgesOnCell, globtolocaledge)# = my_mesh_file["edgesOnCell"][:][:,cells]#, globtolocaledge)
		globalToLocal!(mpasOcean.cellsOnCell, globtolocalcell)# = my_mesh_file["cellsOnCell"][:][:,cells]#, globtolocalcell)
		globalToLocal!(mpasOcean.verticesOnCell, globtolocalvertex)# = my_mesh_file["verticesOnCell"][:][:,cells]#, globtolocalvertex)

		# println("edges on cell post $(mpasOcean.edgesOnCell[1:20])")

		globalToLocal!(mpasOcean.cellsOnEdge, globtolocalcell)# = my_mesh_file["cellsOnEdge"][:][:,edges]#, globtolocalcell)
		ordermapEOE = globalToLocal!(mpasOcean.edgesOnEdge, globtolocaledge)# = my_mesh_file["edgesOnEdge"][:][:,edges], globtolocaledge)
		globalToLocal!(mpasOcean.verticesOnEdge, globtolocalvertex)# = globalToLocal(my_mesh_file["verticesOnEdge"][:][:,edges], globtolocalvertex)

		ordermapCOV = globalToLocal!(mpasOcean.cellsOnVertex, globtolocalcell)# = globalToLocal(my_mesh_file["cellsOnVertex"][:][:,vertices], globtolocalcell)
		globalToLocal!(mpasOcean.edgesOnVertex, globtolocaledge)#
# end
# @timeit tmr "count" begin
		mpasOcean.nEdgesOnCell = dropdims(sum(map(val -> Int64(val != 0), mpasOcean.edgesOnCell), dims=1), dims=1) #base_mesh_file["nEdgesOnCell"][:][cells]
		mpasOcean.nEdgesOnEdge = dropdims(sum(map(val -> Int64(val != 0), mpasOcean.edgesOnEdge), dims=1), dims=1) #base_mesh_file["nEdgesOnCell"][:][cells]
# end
# @timeit tmr "coords" begin

		mpasOcean.xEdge = my_mesh_file["xEdge"][:][edges]
		mpasOcean.yEdge = my_mesh_file["yEdge"][:][edges]


		mpasOcean.latCell = base_mesh_file["latCell"][:][cells]
		mpasOcean.lonCell = base_mesh_file["lonCell"][:][cells]
		mpasOcean.xCell = base_mesh_file["xCell"][:][cells]
		mpasOcean.yCell = base_mesh_file["yCell"][:][cells]
		mpasOcean.areaCell = my_mesh_file["areaCell"][:][cells]

		mpasOcean.dvEdge = my_mesh_file["dvEdge"][:][edges]
		mpasOcean.dcEdge = my_mesh_file["dcEdge"][:][edges]

		mpasOcean.kiteAreasOnVertex = my_mesh_file["kiteAreasOnVertex"][:][:,vertices][ordermapCOV]
		mpasOcean.areaTriangle = my_mesh_file["areaTriangle"][:][vertices]

		mpasOcean.latVertex = base_mesh_file["latVertex"][:][vertices]
		mpasOcean.lonVertex = base_mesh_file["lonVertex"][:][vertices]
		mpasOcean.xVertex = base_mesh_file["xVertex"][:][vertices]
		mpasOcean.yVertex = base_mesh_file["yVertex"][:][vertices]

		# NCDatasets.load!(NCDatasets.variable(my_mesh_file,"weightsOnEdge"),mpasOcean.weightsOnEdge,:,:)[:,edges]
		mpasOcean.weightsOnEdge = my_mesh_file["weightsOnEdge"][:][:,edges][ordermapEOE]


	#         mpasOcean.gridSpacing = mesh_file["gridSpacing"][:]
		mpasOcean.gridSpacing = sqrt.(mpasOcean.areaCell)#mpasOcean.xCell[2] - mpasOcean.xCell[1]

		useGridSpacingMagnitudeDefaultDefinition = true
		if useGridSpacingMagnitudeDefaultDefinition
			mpasOcean.gridSpacingMagnitude = mpasOcean.gridSpacing[0+1]
		else
			mpasOcean.gridSpacingMagnitude = mpasOcean.xCell[1+1] - mpasOcean.xCell[0+1]
		end


		# For a mesh with non-periodic zonal boundaries, adjust the zonal coordinates of the cell centers,
		# vertices and edges.
		dx = mpasOcean.gridSpacingMagnitude # Int64.e. dx = mpasOcean.dcEdge[0]
		dy = sqrt(3.0)/2.0*dx
		if periodicity == "NonPeriodic_x" || periodicity == "NonPeriodic_xy"
			mpasOcean.xCell[:] .-= dx
			mpasOcean.xVertex[:] .-= dx
			mpasOcean.xEdge[:] .-= dx
		end
		if periodicity == "NonPeriodic_y" || periodicity == "NonPeriodic_xy"
			mpasOcean.yCell[:] .-= dy
			mpasOcean.yVertex[:] .-= dy
			mpasOcean.yEdge[:] .-= dy
		end

		# Specify the zonal && meridional extents of the domain.
		mpasOcean.lX = round(maximum(mpasOcean.xCell)) # Int64.e. mpasOcean.lX = sqrt(float(mpasOcean.nCells))*dx
		# Please note that for all of our test problems, the MPAS-O mesh is generated in such a way that
		# mpasOcean.lX Int64.e. the number of cells in the x (|| y) direction times dx has zero fractional part in
		# units of m, for which we can afford to round it to attain perfection.
		mpasOcean.lY = sqrt(3.0)/2.0*mpasOcean.lX # Int64.e. mpasOcean.lY = max(mpasOcean.yVertex)

# end
# @timeit tmr "angleedge" begin
				mpasOcean.angleEdge = my_mesh_file["angleEdge"][:][edges]
				# if true
				mpasOcean.angleEdge[:] = fix_angleEdge(mpasOcean,determineYCellAlongLatitude=false, #true
										   printOutput=false,printRelevantMeshData=false)
				# end

# end
# @timeit tmr "zeros1" begin
	#         mpasOcean.ExactSolutionParameters = DetermineExactSolutionParameters(mpasOcean,false)
		# Define && initialize the following arrays not contained within either the base mesh file || the mesh
		# file.
		mpasOcean.fVertex = zeros(Float64, nVertices)
		mpasOcean.fCell = zeros(Float64, nCells)
		mpasOcean.fEdge = zeros(Float64, nEdges)
		mpasOcean.bottomDepth = zeros(Float64, nCells)
		mpasOcean.bottomDepthEdge = zeros(Float64, nEdges)

		mpasOcean.layerThickness = zeros(Float64, (nCells, mpasOcean.nVertLevels))
		mpasOcean.layerThicknessEdge = zeros(Float64, (nEdges, mpasOcean.nVertLevels))

		DetermineCoriolisParameterAndBottomDepth!(mpasOcean)
# end
# @timeit tmr "max edge z" begin
		mpasOcean.kiteIndexOnCell = zeros(Int64, (nCells,maxEdges))
		mpasOcean.edgeSignOnVertex = zeros(Int8, (nVertices,maxEdges))

		mpasOcean.edgeSignOnCell = zeros(Int8, (nCells,maxEdges))
# end
# @timeit tmr "zeros" begin
		ocn_init_routines_setup_sign_and_index_fields!(mpasOcean)


		mpasOcean.boundaryCell = zeros(Int64, (nCells, mpasOcean.nVertLevels))
		mpasOcean.boundaryEdge = zeros(Int64, (nEdges, mpasOcean.nVertLevels))
		mpasOcean.boundaryVertex = zeros(Int64, (nVertices, mpasOcean.nVertLevels))
		mpasOcean.cellMask = zeros(Int64, (nCells, mpasOcean.nVertLevels))
		mpasOcean.edgeMask = zeros(Int64, (nEdges, mpasOcean.nVertLevels))
		mpasOcean.vertexMask = zeros(Int64, (nVertices, mpasOcean.nVertLevels))

		mpasOcean.maxLevelCell = zeros(Int64, nCells)
		mpasOcean.maxLevelEdgeTop = zeros(Int64, nEdges)
		mpasOcean.maxLevelEdgeBot = zeros(Int64, nEdges)
		mpasOcean.maxLevelVertexTop = zeros(Int64, nVertices)
		mpasOcean.maxLevelVertexBot = zeros(Int64, nVertices)

# end
# @timeit tmr "tendency" begin




		## defining the prognostic variables
		mpasOcean.sshCurrent = zeros(Float64, nCells)
		mpasOcean.sshTendency = zeros(Float64, nCells)

		mpasOcean.normalVelocityCurrent = zeros(Float64, (nEdges, mpasOcean.nVertLevels))
		mpasOcean.normalVelocityTendency = zeros(Float64, (nEdges, mpasOcean.nVertLevels))

		mpasOcean.gravity = 9.8

		# calculate minimum dt based on CFL condition
		mpasOcean.dt = 0.1 * minimum(mpasOcean.dcEdge) / sqrt(mpasOcean.gravity * maximum(mpasOcean.bottomDepth))

		ocn_init_routines_compute_max_level!(mpasOcean)
# end
# show(tmr)
		return mpasOcean
	end
end

#=
function MPAS_Ocean(mesh_directory = "MPAS_O_Shallow_Water/Mesh+Initial_Condition+Registry_Files/Periodic",
					base_mesh_file_name = "base_mesh.nc",
					mesh_file_name = "mesh.nc";
					nvlevels =1,
					periodicity = "Periodic")# where Float64<:AbstractFloat
	return MPAS_Ocean(mesh_directory, base_mesh_file_name, mesh_file_name; nvlevels=nvlevels, periodicity=periodicity, cells=:)
end
function MPAS_Ocean(temp::Bool, mesh_directory = "MPAS_O_Shallow_Water/Mesh+Initial_Condition+Registry_Files/Periodic",
                        base_mesh_file_name = "base_mesh.nc",
                        mesh_file_name = "mesh.nc";
                        periodicity = "Periodic")
    return MPAS_Ocean(mesh_directory, base_mesh_file_name, mesh_file_name, periodicity=periodicity)
end



function MPAS_Ocean(init_state_file_name; periodicity="Periodic")
    mpasOcean = MPAS_Ocean("/", init_state_file_name, init_state_file_name, periodicity=periodicity)
    return mpasOcean
end

=#

#
# function moveArrays!(mpasOcean::MPAS_Ocean, array_type, float_type="default")
#     if float_type=="default"
#         float_type = typeof(mpasOcean.gravity)
#     end
#
#     mainArrayNames = [
#         :sshTendency,
#         :sshCurrent,
#         :nEdgesOnCell,
#         :edgesOnCell,
#         :cellsOnCell,
#         :bottomDepth,
#         :areaCell,
#         :edgeSignOnCell,
#     # ]
#
#     # mainEdgeArrayNames = [
#         :normalVelocityTendency,
#         :normalVelocityCurrent,
#         :cellsOnEdge,
#         :nEdgesOnEdge,
#         :edgesOnEdge,
#         :weightsOnEdge,
#         :fEdge,
#         :dcEdge,
#         :dvEdge,
#         :boundaryEdge,
#         :yEdge,
#         :xEdge,
#         :angleEdge
#     ]
#
#     for (iField, field) in enumerate(mainArrayNames)#fieldnames(typeof(mpasOcean)))
#         # if typeof(mpasOcean).types[iField] == Array
#         if eltype(getfield(mpasOcean, field)) <: AbstractFloat
#             elT = float_type
#         else
#             elT = eltype(getfield(mpasOcean, field))
#         end
#         setfield!(mpasOcean, field, convert(array_type{elT}, getfield(mpasOcean, field)))
#         # end
#     end
#
#
#
# end






function DetermineCoriolisParameterAndBottomDepth!(mpasOcean)
#     alpha0 = mpasOcean.ExactSolutionParameters[0+1]
#     beta0 = mpasOcean.ExactSolutionParameters[2+1]
    f0 = 0.0001 #mpasOcean.ExactSolutionParameters[12+1]
    H = 1000.0 #mpasOcean.ExactSolutionParameters[14+1]
#     if (mpasOcean.myNamelist.config_problem_type == "default"
#         || mpasOcean.myNamelist.config_problem_type == "Barotropic_Tide"
#         || mpasOcean.myNamelist.config_problem_type_Geophysical_Wave)
#         if mpasOcean.myNamelist.config_problem_type == "Planetary_Rossby_Wave"
#             mpasOcean.fCell[:] = f0 + beta0*mpasOcean.yCell[:]
#             mpasOcean.fEdge[:] = f0 + beta0*mpasOcean.yEdge[:]
#             mpasOcean.fVertex[:] = f0 + beta0*mpasOcean.yVertex[:]
#         else
    mpasOcean.fCell[:] .= f0
    mpasOcean.fEdge[:] .= f0
    mpasOcean.fVertex[:] .= f0
#         end
#         if mpasOcean.myNamelist.config_problem_type == "Topographic_Rossby_Wave"
#             mpasOcean.bottomDepth[:] = H + alpha0*yCell[:]
#         else
    mpasOcean.bottomDepth[:] .= H
    for k in 1:mpasOcean.nVertLevels
        mpasOcean.layerThickness[:,k] = mpasOcean.bottomDepth ./ mpasOcean.nVertLevels
    end
    mpasOcean.bottomDepthEdge[:] .= H
    for k in 1:mpasOcean.nVertLevels
        mpasOcean.layerThicknessEdge[:,k] = mpasOcean.bottomDepthEdge ./ mpasOcean.nVertLevels
    end
#         end
#     end
end





function ocn_init_routines_setup_sign_and_index_fields!(mpasOcean)
    for iCell in range(1,mpasOcean.nCells,step=1)
        for j in range(1,mpasOcean.nEdgesOnCell[iCell],step=1)
            iEdge = mpasOcean.edgesOnCell[j,iCell]
            iVertex = mpasOcean.verticesOnCell[j,iCell]
            # Vector points from cell 1 to cell 2
            if mpasOcean.cellsOnEdge[1,iEdge] == iCell
                mpasOcean.edgeSignOnCell[iCell,j] = -1
            else
                mpasOcean.edgeSignOnCell[iCell,j] = 1
            end
            for jj in range(1,mpasOcean.vertexDegree,step=1)
                if mpasOcean.cellsOnVertex[jj,iVertex] == iCell
                    mpasOcean.kiteIndexOnCell[iCell,j] = jj
                end
            end
        end
    end


    for iVertex in range(1,mpasOcean.nVertices,step=1)
        for j in range(1,mpasOcean.vertexDegree,step=1)
            iEdge = mpasOcean.edgesOnVertex[j,iVertex]
            if iEdge != 0
                # Vector points from vertex 1 to vertex 2
                if mpasOcean.verticesOnEdge[1,iEdge] == iVertex
                    mpasOcean.edgeSignOnVertex[iVertex,j] = -1
                else
                    mpasOcean.edgeSignOnVertex[iVertex,j] = 1
                end
            end
        end
    end
end



function ocn_init_routines_compute_max_level!(mpasOcean)
    mpasOcean.maxLevelCell[:] .= mpasOcean.nVertLevels-1

    for iEdge in 1:mpasOcean.nEdges
        iCell1 = mpasOcean.cellsOnEdge[1,iEdge]
        iCell2 = mpasOcean.cellsOnEdge[2,iEdge]
        if iCell1 == 0 || iCell2 == 0
            mpasOcean.boundaryEdge[iEdge,:] .= 1.0
            mpasOcean.maxLevelEdgeTop[iEdge] = -1
            if iCell1 == 0
                mpasOcean.maxLevelEdgeBot[iEdge] = mpasOcean.maxLevelCell[iCell2]
            elseif iCell2 == 0
                mpasOcean.maxLevelEdgeBot[iEdge] = mpasOcean.maxLevelCell[iCell1]
            end
        elseif ! ( iCell1 == 0 || iCell2 == 0 )
            mpasOcean.maxLevelEdgeTop[iEdge] = min(mpasOcean.maxLevelCell[iCell1], mpasOcean.maxLevelCell[iCell2])
            mpasOcean.maxLevelEdgeBot[iEdge] = max(mpasOcean.maxLevelCell[iCell1], mpasOcean.maxLevelCell[iCell2])
        end
    end

    for iVertex in 1:mpasOcean.nVertices
        iCell1 = mpasOcean.cellsOnVertex[1,iVertex]
        if iCell1 == 0
            mpasOcean.maxLevelVertexBot[iVertex] = -1
            mpasOcean.maxLevelVertexTop[iVertex] = -1
        else
            mpasOcean.maxLevelVertexBot[iVertex] = mpasOcean.maxLevelCell[iCell1]
            mpasOcean.maxLevelVertexTop[iVertex] = mpasOcean.maxLevelCell[iCell1]
        end

        for Int64 in 1:mpasOcean.vertexDegree
            iCell = mpasOcean.cellsOnVertex[Int64, iVertex]
            if iCell == 0
                mpasOcean.maxLevelVertexBot[iVertex] = max(mpasOcean.maxLevelVertexBot[iVertex], -1)
                mpasOcean.maxLevelVertexBot[iVertex] = min(mpasOcean.maxLevelVertexTop[iVertex], -1)
            else
                mpasOcean.maxLevelVertexBot[iVertex] = max(mpasOcean.maxLevelVertexBot[iVertex],
                                                            mpasOcean.maxLevelCell[iCell])
                mpasOcean.maxLevelVertexTop[iVertex] = min(mpasOcean.maxLevelVertexTop[iVertex],
                                                            mpasOcean.maxLevelCell[iCell])
            end
        end
    end

    # determine_boundaryEdge_Generalized_Method = true
	#
    # if determine_boundaryEdge_Generalized_Method
    #     mpasOcean.boundaryEdge[:,:] .= 1
    # end
    mpasOcean.edgeMask[:,:] .= 0

    for iEdge in 1:mpasOcean.nEdges
        index_UpperLimit = mpasOcean.maxLevelEdgeTop[iEdge]
        if index_UpperLimit > -1
#             if determine_boundaryEdge_Generalized_Method
#                 mpasOcean.boundaryEdge[iEdge,1:index_UpperLimit+1] .= 0
#             end
            mpasOcean.edgeMask[iEdge,1:index_UpperLimit+1] .= 1
        end
    end

    for iEdge in 1:mpasOcean.nEdges
        iCell1 = mpasOcean.cellsOnEdge[1, iEdge]
        iCell2 = mpasOcean.cellsOnEdge[2, iEdge]

        iVertex1 = mpasOcean.verticesOnEdge[1, iEdge]
        iVertex2 = mpasOcean.verticesOnEdge[2, iEdge]

        if mpasOcean.boundaryEdge[iEdge] == 1
            if iCell1 != 0
                mpasOcean.boundaryCell[iCell1] = 1
            end
            if iCell2 != 0
                mpasOcean.boundaryCell[iCell2] = 1
            end
            mpasOcean.boundaryVertex[iVertex1] = 1
            mpasOcean.boundaryVertex[iVertex2] = 1
        end
    end


    for iCell in 1:mpasOcean.nCells
        for k in 1:mpasOcean.nVertLevels
            if mpasOcean.maxLevelCell[iCell] >= k
                mpasOcean.cellMask[iCell, k] = 1
            end
        end
    end
    for iVertex in 1:mpasOcean.nVertices
        for k in 1:mpasOcean.nVertLevels
            if mpasOcean.maxLevelVertexBot[iVertex] >= k
                mpasOcean.vertexMask[iVertex, k] = 1
            end
        end
    end

    mpasOcean.nNonPeriodicBoundaryEdges = 0.0
    for iEdge in 1:mpasOcean.nEdges
        if mpasOcean.boundaryEdge[iEdge, 1] == 1.0
            mpasOcean.nNonPeriodicBoundaryEdges += 1
        end
    end

	mpasOcean.nNonPeriodicBoundaryVertices = 0.0
    for iVertex in 1:mpasOcean.nVertices
        if mpasOcean.boundaryVertex[iVertex, 1] == 1.0
            mpasOcean.nNonPeriodicBoundaryVertices += 1
        end
    end

	mpasOcean.nNonPeriodicBoundaryCells = 0.0
    for iCell in 1:mpasOcean.nCells
        if mpasOcean.boundaryCell[iCell, 1] == 1.0
            mpasOcean.nNonPeriodicBoundaryCells += 1
        end
    end
end
