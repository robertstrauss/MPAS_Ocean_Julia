import NCDatasets
import CUDA
# import KernelAbstractions

# include("Namelist.jl")
include("fixAngleEdge.jl")
include("BoundaryCondition.jl")




mutable struct MPAS_Ocean{F<:AbstractFloat}
    
#     boundaryConditions::Array{BoundaryCondition}
    
    # prognostic variables
    sshCurrent::AbstractArray{F}   #where F <: AbstractFloat # sea surface height (cell-centered)
    sshTendency::AbstractArray{F}   #where F <: AbstractFloat # tendency (cell-centered)
    
    normalVelocityCurrent::AbstractArray{F}   #where F <: AbstractFloat # group velocity normal to mesh edges (edge-centered)
    normalVelocityTendency::AbstractArray{F}   #where F <: AbstractFloat # tendency (edge-centered)
    
    
    bottomDepth::AbstractArray{F}   #where F <: AbstractFloat # bathymetry (cell-centered)
    gravity::F
    
    dt::F # timestep
    
    
    
    ### mesh information
    
    
    ## cell-centered arrays
    nCells::Integer # number of cells
    cellsOnCell::AbstractArray{I} where I <: Integer # indices of neighbor cells (nCells, nEdgesOnCell[i])
    edgesOnCell::AbstractArray{I} where I <: Integer # indices of edges of cell (nCells, nEdgesOnCell[i])
    verticesOnCell::AbstractArray{I} where I <: Integer # indices of vertices of cell (nCells, nEdgesOnCell[i])
    kiteIndexOnCell::AbstractArray{I} where I <: Integer # index of kite shape on cell connecting vertex, midpoint of edge, center of cell, and midpoint of other edge (nCells, nEdgesOnCell[i])
    nEdgesOnCell::AbstractArray{I} where I <: Integer # number of edges a cell has (cell-centered)
    edgeSignOnCell::AbstractArray{Int8} # orientation of edge relative to cell (nCells, nEdgesOnCell[i])
    latCell::AbstractArray{F}   #where F <: AbstractFloat # latitude of cell (cell-centered)
    lonCell::AbstractArray{F}   #where F <: AbstractFloat # longitude of cell (cell-centered)
    xCell::AbstractArray{F}   #where F <: AbstractFloat # x coordinate of cell (cell-centered)
    yCell::AbstractArray{F}   #where F <: AbstractFloat # y coordinate of cell (cell-centered)
    areaCell::AbstractArray{F}   #where F <: AbstractFloat # area of cell (cell-centered)
    fCell::AbstractArray{F}   #where F <: AbstractFloat # coriolis parameter (cell-centered)
    maxLevelCell::AbstractArray{I} where I <: Integer
    gridSpacing::AbstractArray{F}   #where F <: AbstractFloat
    boundaryCell::AbstractArray{I} where I <: Integer # 0 for inner cells, 1 for boundary cells (cell-centered)
    
    
    ## edge-centered arrays
    nEdges::Integer
    cellsOnEdge::AbstractArray{I} where I <: Integer # indices of the two cells on this edge (nEdges, 2)
    edgesOnEdge::AbstractArray{I} where I <: Integer
    verticesOnEdge::AbstractArray{I} where I <: Integer
    nEdgesOnEdge::AbstractArray{I} where I <: Integer
    xEdge::AbstractArray{F}   #where F <: AbstractFloat
    yEdge::AbstractArray{F}   #where F <: AbstractFloat
    dvEdge::AbstractArray{F}   #where F <: AbstractFloat
    dcEdge::AbstractArray{F}   #where F <: AbstractFloat
    fEdge::AbstractArray{F}   #where F <: AbstractFloat # coriolis parameter
    angleEdge::AbstractArray{F}   #where F <: AbstractFloat
    weightsOnEdge::AbstractArray{F}   #where F <: AbstractFloat # coeffecients of norm vels of surrounding edges in linear combination to compute tangential velocity
    maxLevelEdgeTop::AbstractArray{I} where I <: Integer
    maxLevelEdgeBot::AbstractArray{I} where I <: Integer
    boundaryEdge::AbstractArray{I} where I <: Integer
    edgeMask::AbstractArray{I} where I <: Integer
    
    
    ## vertex-centered arrays
    nVertices::Integer
    latVertex::AbstractArray{F}   #where F <: AbstractFloat
    lonVertex::AbstractArray{F}   #where F <: AbstractFloat
    xVertex::AbstractArray{F}   #where F <: AbstractFloat
    yVertex::AbstractArray{F}   #where F <: AbstractFloat
    vertexDegree::Integer
    cellsOnVertex::AbstractArray{I} where I <: Integer
    edgesOnVertex::AbstractArray{I} where I <: Integer
    edgeSignOnVertex::AbstractArray{I} where I <: Integer
    fVertex::AbstractArray{F}   #where F <: AbstractFloat
    areaTriangle::AbstractArray{F}   #where F <: AbstractFloat
    kiteAreasOnVertex::AbstractArray{F}   #where F <: AbstractFloat
    maxLevelVertexTop::AbstractArray{I} where I <: Integer
    maxLevelVertexBot::AbstractArray{I} where I <: Integer
    boundaryVertex::AbstractArray{I} where I <: Integer

    
    
    nVertLevels::Integer


#     ExactSolutionParameters::AbstractArray{F}   #where F <: AbstractFloat
#     myNamelist::Namelist

    gridSpacingMagnitude::F

    # maximum values of x and y (size of mesh)
    lX::F
    lY::F

    nNonPeriodicBoundaryCells::Integer
    nNonPeriodicBoundaryEdges::Integer
    nNonPeriodicBoundaryVertices::Integer




    
    # load state from file
#     function MPAS_Ocean(mesh_directory::String, base_mesh_file_name, mesh_file_name; periodicity::String)
#         return MPAS_Ocean(mesh_directory=mesh_directory, base_mesh_file_name=base_mesh_file_name, mesh_file_name=mesh_file_name, periodicity=periodicity)
#     end
    function MPAS_Ocean(mesh_directory = "MPAS_O_Shallow_Water/Mesh+Initial_Condition+Registry_Files/Periodic",
                        base_mesh_file_name = "base_mesh.nc",
                        mesh_file_name = "mesh.nc"; periodicity = "Periodic")
        return MPAS_Ocean{Float64}( mesh_directory,
                        base_mesh_file_name,
                        mesh_file_name;
                        periodicity = "Periodic")
    end
    function MPAS_Ocean{F}(mesh_directory = "MPAS_O_Shallow_Water/Mesh+Initial_Condition+Registry_Files/Periodic",
                        base_mesh_file_name = "base_mesh.nc",
                        mesh_file_name = "mesh.nc";
                        periodicity = "Periodic") where F<:AbstractFloat
        mpasOcean = new{F}()
        
        


#        mpasOcean.boundaryConditions = [BoundaryCondition()]
        
        
        ### Mesh stuff
#         mpasOcean.myNamelist = Namelist()
        
        base_mesh_file = NCDatasets.Dataset("$mesh_directory/$base_mesh_file_name", "r", format=:netcdf4_classic)
        mesh_file = NCDatasets.Dataset("$mesh_directory/$mesh_file_name", "r", format=:netcdf4_classic)

        maxEdges = mesh_file.dim["maxEdges"]

        ## defining the mesh
        mpasOcean.nCells = mesh_file.dim["nCells"]
        mpasOcean.nEdges = mesh_file.dim["nEdges"]
        mpasOcean.nVertices = mesh_file.dim["nVertices"]
        mpasOcean.vertexDegree = mesh_file.dim["vertexDegree"]
        mpasOcean.nVertLevels = 1

        # will be incremented in ocn_init_routines_compute_max_level
        mpasOcean.nNonPeriodicBoundaryEdges = 0
        mpasOcean.nNonPeriodicBoundaryVertices = 0
        mpasOcean.nNonPeriodicBoundaryCells = 0

        nCells = mpasOcean.nCells
        nEdges = mpasOcean.nEdges
        nVertices = mpasOcean.nVertices

        # Choose my_mesh_file_name to be either base_mesh_file_name or mesh_file_name
        my_mesh_file_name = mesh_file_name
        # Choose my_mesh_file to be either base_mesh_file or mesh_file
        my_mesh_file = mesh_file

        mpasOcean.cellsOnCell = my_mesh_file["cellsOnCell"][:]
        mpasOcean.nEdgesOnCell = base_mesh_file["nEdgesOnCell"][:]
        mpasOcean.edgesOnCell = my_mesh_file["edgesOnCell"][:]
        mpasOcean.verticesOnCell = my_mesh_file["verticesOnCell"][:]

        mpasOcean.cellsOnEdge = my_mesh_file["cellsOnEdge"][:]
        mpasOcean.nEdgesOnEdge = base_mesh_file["nEdgesOnEdge"][:]
        mpasOcean.edgesOnEdge = my_mesh_file["edgesOnEdge"][:]
        mpasOcean.verticesOnEdge = my_mesh_file["verticesOnEdge"][:]

        mpasOcean.cellsOnVertex = my_mesh_file["cellsOnVertex"][:]
        mpasOcean.edgesOnVertex = my_mesh_file["edgesOnVertex"][:]
        #


        mpasOcean.angleEdge = my_mesh_file["angleEdge"][:]
        # if true
        mpasOcean.angleEdge[:] = fix_angleEdge(mesh_directory,my_mesh_file_name,determineYCellAlongLatitude=true,
                                   printOutput=false,printRelevantMeshData=false)
        # end



        mpasOcean.xEdge = my_mesh_file["xEdge"][:]
        mpasOcean.yEdge = my_mesh_file["yEdge"][:]


        mpasOcean.latCell = base_mesh_file["latCell"][:]
        mpasOcean.lonCell = base_mesh_file["lonCell"][:]
        mpasOcean.xCell = base_mesh_file["xCell"][:]
        mpasOcean.yCell = base_mesh_file["yCell"][:]

        mpasOcean.areaCell = my_mesh_file["areaCell"][:]
        mpasOcean.dvEdge = my_mesh_file["dvEdge"][:]
        mpasOcean.dcEdge = my_mesh_file["dcEdge"][:]

        mpasOcean.kiteAreasOnVertex = my_mesh_file["kiteAreasOnVertex"][:]
        mpasOcean.areaTriangle = my_mesh_file["areaTriangle"][:]


        mpasOcean.latVertex = base_mesh_file["latVertex"][:]
        mpasOcean.lonVertex = base_mesh_file["lonVertex"][:]
        mpasOcean.xVertex = base_mesh_file["xVertex"][:]
        mpasOcean.yVertex = base_mesh_file["yVertex"][:]


        mpasOcean.gridSpacing = mesh_file["gridSpacing"][:]

        useGridSpacingMagnitudeDefaultDefinition = true
        if useGridSpacingMagnitudeDefaultDefinition
            mpasOcean.gridSpacingMagnitude = mpasOcean.gridSpacing[0+1]
        else
            mpasOcean.gridSpacingMagnitude = mpasOcean.xCell[1+1] - mpasOcean.xCell[0+1]
        end


        # For a mesh with non-periodic zonal boundaries, adjust the zonal coordinates of the cell centers,
        # vertices and edges.
        dx = mpasOcean.gridSpacingMagnitude # i.e. dx = mpasOcean.dcEdge[0]
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
        mpasOcean.lX = round(maximum(mpasOcean.xCell)) # i.e. mpasOcean.lX = sqrt(float(mpasOcean.nCells))*dx
        # Please note that for all of our test problems, the MPAS-O mesh is generated in such a way that
        # mpasOcean.lX i.e. the number of cells in the x (|| y) direction times dx has zero fractional part in
        # units of m, for which we can afford to round it to attain perfection.
        mpasOcean.lY = sqrt(3.0)/2.0*mpasOcean.lX # i.e. mpasOcean.lY = max(mpasOcean.yVertex)


#         mpasOcean.ExactSolutionParameters = DetermineExactSolutionParameters(mpasOcean,false)
        # Define && initialize the following arrays not contained within either the base mesh file || the mesh
        # file.
        mpasOcean.fVertex = zeros(F, nVertices)
        mpasOcean.fCell = zeros(F, nCells)
        mpasOcean.fEdge = zeros(F, nEdges)
        mpasOcean.bottomDepth = zeros(F, nCells)
        DetermineCoriolisParameterAndBottomDepth!(mpasOcean)

        mpasOcean.kiteIndexOnCell = zeros(Int64, (nCells,maxEdges))
        mpasOcean.edgeSignOnVertex = zeros(Int8, (nVertices,maxEdges))


        mpasOcean.weightsOnEdge = my_mesh_file["weightsOnEdge"][:]

        mpasOcean.edgeSignOnCell = zeros(Int8, (nCells,maxEdges))
        ocn_init_routines_setup_sign_and_index_fields!(mpasOcean)


        mpasOcean.boundaryCell = zeros(Int8, (nCells, mpasOcean.nVertLevels))
        mpasOcean.boundaryEdge = zeros(Int8, (nEdges, mpasOcean.nVertLevels))
        mpasOcean.boundaryVertex = zeros(Int8, (nVertices, mpasOcean.nVertLevels))
        mpasOcean.edgeMask = zeros(Int8, nEdges)

        mpasOcean.maxLevelCell = zeros(Int64, nCells)
        mpasOcean.maxLevelEdgeTop = zeros(Int64, nEdges)
        mpasOcean.maxLevelEdgeBot = zeros(Int64, nEdges)
        mpasOcean.maxLevelVertexTop = zeros(Int64, nVertices)
        mpasOcean.maxLevelVertexBot = zeros(Int64, nVertices)


        
        
        
        
        ## defining the prognostic variables
        mpasOcean.sshCurrent = zeros(F, nCells)
        mpasOcean.sshTendency = zeros(F, nCells)

        mpasOcean.normalVelocityCurrent = zeros(F, (nEdges))
        mpasOcean.normalVelocityTendency = zeros(F, (nEdges))
        
        mpasOcean.gravity = 9.8
        
        # calculate minimum dt based on CFL condition
        mpasOcean.dt = 0.1 * minimum(mpasOcean.dcEdge) / sqrt(mpasOcean.gravity * maximum(mpasOcean.bottomDepth))
        
        ocn_init_routines_compute_max_level!(mpasOcean)

        return mpasOcean
    end
end
function MPAS_Ocean(temp::Bool, mesh_directory = "MPAS_O_Shallow_Water/Mesh+Initial_Condition+Registry_Files/Periodic",
                        base_mesh_file_name = "base_mesh.nc",
                        mesh_file_name = "mesh.nc";
                        periodicity = "Periodic")
    return MPAS_Ocean(mesh_directory, base_mesh_file_name, mesh_file_name, periodicity=periodicity)
end








function moveArrays!(mpasOcean::MPAS_Ocean, array_type, float_type="default")
    if float_type=="default"
        float_type = typeof(mpasOcean.gravity)
    end
    
    mainArrayNames = [
        :sshTendency,
        :sshCurrent,
        :nEdgesOnCell,
        :edgesOnCell,
        :cellsOnCell,
        :bottomDepth,
        :areaCell,
        :edgeSignOnCell,
    # ]
    
    # mainEdgeArrayNames = [
        :normalVelocityTendency,
        :normalVelocityCurrent,
        :cellsOnEdge,
        :nEdgesOnEdge,
        :edgesOnEdge,
        :weightsOnEdge,
        :fEdge,
        :dcEdge,
        :dvEdge,
        :boundaryEdge,
        :yEdge,
        :xEdge,
        :angleEdge
    ]
    
    for (iField, field) in enumerate(mainArrayNames)#fieldnames(typeof(mpasOcean)))
        # if typeof(mpasOcean).types[iField] == AbstractArray
        if eltype(getfield(mpasOcean, field)) <: AbstractFloat
            elT = float_type
        else
            elT = eltype(getfield(mpasOcean, field))
        end
        setfield!(mpasOcean, field, convert(array_type{elT}, getfield(mpasOcean, field)))
        # end
    end
    

    
end






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
#         end
#     end
end





function ocn_init_routines_setup_sign_and_index_fields!(mpasOcean)
    for iCell in range(1,mpasOcean.nCells,step=1)
        for i in range(1,mpasOcean.nEdgesOnCell[iCell],step=1)
            iEdge = mpasOcean.edgesOnCell[i,iCell]
            iVertex = mpasOcean.verticesOnCell[i,iCell]
            # Vector points from cell 1 to cell 2
            if mpasOcean.cellsOnEdge[1,iEdge] == iCell
                mpasOcean.edgeSignOnCell[iCell,i] = -1
            else
                mpasOcean.edgeSignOnCell[iCell,i] = 1
            end
            for j in range(1,mpasOcean.vertexDegree,step=1)
                if mpasOcean.cellsOnVertex[j,iVertex] == iCell
                    mpasOcean.kiteIndexOnCell[iCell,i] = j
                end
            end
        end
    end
#     printKiteIndexOnCell = false
#     if printKiteIndexOnCell
#         for iCell in range(1,mpasOcean.nCells,step=1)
#             for i in range(1,mpasOcean.nEdgesOnCell[iCell],step=1)
#                 println("On edge $i of cell $iCell, kiteIndexOnCell is $(mpasOcean.kiteIndexOnCell[iCell,i]).")
#             end
#         end
#     end

    for iVertex in range(1,mpasOcean.nVertices,step=1)
        for i in range(1,mpasOcean.vertexDegree,step=1)
            iEdge = mpasOcean.edgesOnVertex[i,iVertex]
            if iEdge != 0
                # Vector points from vertex 1 to vertex 2
                if mpasOcean.verticesOnEdge[1,iEdge] == iVertex
                    mpasOcean.edgeSignOnVertex[iVertex,i] = -1
                else
                    mpasOcean.edgeSignOnVertex[iVertex,i] = 1
                end
            end
        end
    end
end



function ocn_init_routines_compute_max_level!(mpasOcean)
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

        for i in 1:mpasOcean.vertexDegree
            iCell = mpasOcean.cellsOnVertex[i, iVertex]
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

    determine_boundaryEdge_Generalized_Method = true

    if determine_boundaryEdge_Generalized_Method
        mpasOcean.boundaryEdge[:,:] .= 1
    end
    mpasOcean.edgeMask[:,:] .= 0

    for iEdge in 1:mpasOcean.nEdges
        index_UpperLimit = mpasOcean.maxLevelEdgeTop[iEdge]
        if index_UpperLimit > -1
            if determine_boundaryEdge_Generalized_Method
                mpasOcean.boundaryEdge[iEdge,1:index_UpperLimit+1] .= 0
            end
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


    for iEdge in 1:mpasOcean.nEdges
        if mpasOcean.boundaryEdge[iEdge, 1] == 1.0
            mpasOcean.nNonPeriodicBoundaryEdges += 1
        end
    end
    for iVertex in 1:mpasOcean.nVertices
        if mpasOcean.boundaryVertex[iVertex, 1] == 1.0
            mpasOcean.nNonPeriodicBoundaryVertices += 1
        end
    end
    for iCell in 1:mpasOcean.nCells
        if mpasOcean.boundaryCell[iCell, 1] == 1.0
            mpasOcean.nNonPeriodicBoundaryCells += 1
        end
    end
end


# function DetermineExactSolutionParameters(mpasOcean,printPhaseSpeedOfWaveModes)
#     angleEdge = 0.0
#     beta0 = mpasOcean.myNamelist.config_meridional_gradient_of_Coriolis_parameter
#     c = mpasOcean.myNamelist.config_phase_speed_of_coastal_Kelvin_wave
#     dx = mpasOcean.gridSpacingMagnitude # i.e. dx = mpasOcean.dcEdge[0+1]
#     dy = sqrt(3.0)/2.0*dx
#     etaHat1 = mpasOcean.myNamelist.config_surface_elevation_amplitude
#     if mpasOcean.myNamelist.config_problem_type == "Viscous_Burgers_Equation"
#         etaHat2 = 0.0
#     else
#         etaHat2 = 2.0*etaHat1
#     end
#     f0 = mpasOcean.myNamelist.config_Coriolis_parameter
#     g = mpasOcean.myNamelist.config_gravity
#     H = mpasOcean.myNamelist.config_mean_depth
#     if mpasOcean.myNamelist.config_problem_type == "Topographic_Rossby_Wave"
#         alpha0 = beta0*H/f0
#         # In the Northern Hemisphere where f0 > 0, the topographic Rossby wave travels with the shallower water on
#         # its right. Hence if alpha0 > 0 i.e. the ocean depth increases northward, the topographic Rossby wave
#         # will propagate eastward else it will propagate westward.
#         mpasOcean.myNamelist.config_bottom_slope = alpha0
#     else
#         alpha0 = mpasOcean.myNamelist.config_bottom_slope # alpha0 = 0.0
#     end
#     if mpasOcean.myNamelist.config_problem_type == "Barotropic_Tide"
#         kX1 = 2.5*pi/mpasOcean.lX
#         kY1 = 0.0
#         kX2 = 4.5*pi/mpasOcean.lX
#         kY2 = 0.0
#     elseif mpasOcean.myNamelist.config_problem_type == "Viscous_Burgers_Equation"
#         kX1 = 0.0
#         kY1 = 0.0
#         kX2 = 0.0
#         kY2 = 0.0
#     else
#         kX1 = 2.0*pi/mpasOcean.lX
#         kY1 = 2.0*pi/mpasOcean.lY
#         kX2 = 2.0*kX1
#         kY2 = 2.0*kY1
#     end
#     lX = mpasOcean.lX
#     lY = mpasOcean.lY
#     R = mpasOcean.myNamelist.config_radius_of_deformation
#     Req = mpasOcean.myNamelist.config_equatorial_radius_of_deformation
#     LengthScale = sqrt(c/beta0)
#     TimeScale = 1.0/sqrt(beta0*c)
#     VelocityScale = c
#     SurfaceElevationScale = c^2.0/g
#     if (mpasOcean.myNamelist.config_problem_type == "default"
#         || mpasOcean.myNamelist.config_problem_type == "Coastal_Kelvin_Wave")
#         omega1 = -c*kY1
#         omega2 = -c*kY2
#     elseif mpasOcean.myNamelist.config_problem_type == "Inertia_Gravity_Wave"
#         omega1 = sqrt(g*H*(kX1^2.0 + kY1^2.0) + f0^2.0)
#         omega2 = sqrt(g*H*(kX2^2.0 + kY2^2.0) + f0^2.0)
#     elseif mpasOcean.myNamelist.config_problem_type == "Planetary_Rossby_Wave"
#         omega1 = -beta0*R^2.0*kX1/(1.0 + R^2.0*(kX1^2.0 + kY1^2.0))
#         omega2 = -beta0*R^2.0*kX2/(1.0 + R^2.0*(kX2^2.0 + kY2^2.0))
#     elseif mpasOcean.myNamelist.config_problem_type == "Topographic_Rossby_Wave"
#         omega1 = alpha0*g/f0*kX1/(1.0 + R^2.0*(kX1^2.0 + kY1^2.0))
#         omega2 = alpha0*g/f0*kX2/(1.0 + R^2.0*(kX2^2.0 + kY2^2.0))
#     elseif mpasOcean.myNamelist.config_problem_type == "Equatorial_Kelvin_Wave"
#         omega1 = c*kX1
#         omega2 = c*kX2
#     elseif mpasOcean.myNamelist.config_problem_type == "Equatorial_Yanai_Wave"
#         omega1 = GWESST.DetermineEquatorialYanaiWaveNonDimensionalAngularFrequency(LengthScale*kX1)/TimeScale
#         omega2 = GWESST.DetermineEquatorialYanaiWaveNonDimensionalAngularFrequency(LengthScale*kX2)/TimeScale
#     elseif mpasOcean.myNamelist.config_problem_type == "Equatorial_Rossby_Wave"
#         omega1 = GWESST.DetermineEquatorialRossbyAndInertiaGravityWaveNonDimensionalAngularFrequency(
#                  "Equatorial_Rossby_Wave",LengthScale*kX1,m=1)/TimeScale
#         omega2 = GWESST.DetermineEquatorialRossbyAndInertiaGravityWaveNonDimensionalAngularFrequency(
#                  "Equatorial_Rossby_Wave",LengthScale*kX2,m=1)/TimeScale
#     elseif mpasOcean.myNamelist.config_problem_type == "Equatorial_Inertia_Gravity_Wave"
#         omega1 = GWESST.DetermineEquatorialRossbyAndInertiaGravityWaveNonDimensionalAngularFrequency(
#                  "Equatorial_Inertia_Gravity_Wave",LengthScale*kX1,m=2)/TimeScale
#         omega2 = GWESST.DetermineEquatorialRossbyAndInertiaGravityWaveNonDimensionalAngularFrequency(
#                  "Equatorial_Inertia_Gravity_Wave",LengthScale*kX2,m=2)/TimeScale
#     elseif mpasOcean.myNamelist.config_problem_type == "Barotropic_Tide"
#         omega1 = sqrt(g*H*kX1^2.0 + f0^2.0)
#         omega2 = sqrt(g*H*kX2^2.0 + f0^2.0)
#     elseif (mpasOcean.myNamelist.config_problem_type == "Diffusion_Equation"
#           || mpasOcean.myNamelist.config_problem_type == "Viscous_Burgers_Equation")
#         omega1 = 0.0
#         omega2 = 0.0
#     end
#     if (mpasOcean.myNamelist.config_problem_type == "default"
#         || mpasOcean.myNamelist.config_problem_type == "Coastal_Kelvin_Wave")
#         cX1 = 0.0
#         cY1 = -c # cY1 = omega1/kY1 = -c
#         cX2 = 0.0
#         cY2 = -c # cY2 = omega2/kY2 = -c
#     elseif mpasOcean.myNamelist.config_problem_type_Equatorial_Wave
#         if mpasOcean.myNamelist.config_problem_type == "Equatorial_Kelvin_Wave"
#             cX1 = c
#             cY1 = 0.0
#             cX2 = c
#             cY2 = 0.0
#         else
#             cX1 = omega1/kX1
#             cY1 = 0.0
#             cX2 = omega2/kX2
#             cY2 = 0.0
#         end
#     elseif mpasOcean.myNamelist.config_problem_type == "Barotropic_Tide"
#         cX1 = omega1/kX1
#         cY1 = 0.0
#         cX2 = omega2/kX2
#         cY2 = 0.0
#     elseif (mpasOcean.myNamelist.config_problem_type == "Diffusion_Equation"
#           || mpasOcean.myNamelist.config_problem_type == "Viscous_Burgers_Equation")
#         cX1 = 0.0
#         cY1 = 0.0
#         cX2 = 0.0
#         cY2 = 0.0
#     else
#         cX1 = omega1/kX1
#         cY1 = omega1/kY1
#         cX2 = omega2/kX2
#         cY2 = omega2/kY2
#     end
#     if ((mpasOcean.myNamelist.config_problem_type_Geophysical_Wave
#          || mpasOcean.myNamelist.config_problem_type == "Barotropic_Tide") && printPhaseSpeedOfWaveModes)
#         println("The zonal component of the phase speed of the first wave mode is %.4g." %cX1)
#         println("The meridional component of the phase speed of the first wave mode is %.4g." %cY1)
#         println("The zonal component of the phase speed of the second wave mode is %.4g." %cX2)
#         println("The meridional component of the phase speed of the second wave mode is %.4g." %cY2)
#     end
#     kappaX = mpasOcean.myNamelist.config_zonal_diffusivity
#     kappaY = mpasOcean.myNamelist.config_meridional_diffusivity
#     kappa1 = kappaX*kX1^2.0 + kappaY*kY1^2.0
#     kappa2 = kappaX*kX2^2.0 + kappaY*kY2^2.0
#     uL = mpasOcean.myNamelist.config_viscous_Burgers_zonal_velocity_left
#     uR = mpasOcean.myNamelist.config_viscous_Burgers_zonal_velocity_right
#     s = mpasOcean.myNamelist.config_viscous_Burgers_shock_speed
#     mpasOcean.myNamelist.config_viscous_Burgers_zonal_offset = 0.25*lX
#     x0 = mpasOcean.myNamelist.config_viscous_Burgers_zonal_offset
#     ExactSolutionParameters = [alpha0,angleEdge,beta0,c,cX1,cX2,cY1,cY2,dx,dy,etaHat1,etaHat2,f0,g,H,kX1,kX2,kY1,
#                                kY2,lX,lY,omega1,omega2,R,Req,LengthScale,TimeScale,VelocityScale,
#                                SurfaceElevationScale,kappaX,kappaY,kappa1,kappa2,uL,uR,s,x0]
#     return ExactSolutionParameters
# end



