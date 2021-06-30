import NCDatasets
using KernelAbstractions

include("Namelist.jl")
include("fixAngleEdge.jl")
gravity = 9.8


# Structure for storing state of ocean
mutable struct MPAS_Ocean
    nCells::Integer # number of cells
    # all below have length nCells, index i in any list is about cell #i
    sshCurrent::AbstractArray{F} where F <: AbstractFloat # 1d list of sea surface heights of cells
    sshNew::AbstractArray{F} where F <: AbstractFloat # 1d future ssh
    bottomDepth::AbstractArray{F} where F <: AbstractFloat # 1d list of bottomDepth/ocean floor topography of cells
    cellsOnCell::AbstractArray{I} where I <: Integer # 2d (nCells x nEdgesOnCell[i]) indices of the cells neighboring this one
    edgesOnCell::AbstractArray{I} where I <: Integer # 2d (nCells x nEdgesOnCell[i]) indices of the edges/normal velocities bordering this one
    verticesOnCell::AbstractArray{I} where I <: Integer # 2d (nCells x nEdgesOnCell[i]) indices of the edges/normal velocities bordering this one
    kiteIndexOnCell::AbstractArray{I} where I <: Integer
    nEdgesOnCell::AbstractArray{I} where I <: Integer # 1d  list of the number of edges a cell has
    edgeSignOnCell::AbstractArray{Int8} # 1d list of sign of norm vel
    latCell::AbstractArray{F} where F <: AbstractFloat
    lonCell::AbstractArray{F} where F <: AbstractFloat
    xCell::AbstractArray{F} where F <: AbstractFloat
    yCell::AbstractArray{F} where F <: AbstractFloat
    areaCell::AbstractArray{F} where F <: AbstractFloat
    fCell::AbstractArray{F} where F <: AbstractFloat
    maxLevelCell::AbstractArray{I} where I <: Integer
    gridSpacing::AbstractArray{F} where F <: AbstractFloat
    boundaryCell::AbstractArray{I} where I <: Integer

    nEdges::Integer
    # all below have length nEdges, && index i in any list is about edge #i
    normalVelocityCurrent::AbstractArray{F} where F <: AbstractFloat # 1d list of norm vel
    normalVelocityNew::AbstractArray{F} where F <: AbstractFloat # 1d list of future norm vel
    cellsOnEdge::AbstractArray{I} where I <: Integer # 2d (nEdges x 2) the indexes of the two cells forming this edge
    edgesOnEdge::AbstractArray{I} where I <: Integer
    verticesOnEdge::AbstractArray{I} where I <: Integer
    nEdgesOnEdge::AbstractArray{I} where I <: Integer
    xEdge::AbstractArray{F} where F <: AbstractFloat
    yEdge::AbstractArray{F} where F <: AbstractFloat
    dvEdge::AbstractArray{F} where F <: AbstractFloat
    dcEdge::AbstractArray{F} where F <: AbstractFloat
    fEdge::AbstractArray{F} where F <: AbstractFloat # coriolis parameter
    angleEdge::AbstractArray{F} where F <: AbstractFloat
    weightsOnEdge::AbstractArray{F} where F <: AbstractFloat # weighting for computing tangential velocity
    maxLevelEdgeTop::AbstractArray{I} where I <: Integer
    maxLevelEdgeBot::AbstractArray{I} where I <: Integer
    boundaryEdge::AbstractArray{I} where I <: Integer
    edgeMask::AbstractArray{I} where I <: Integer


    nVertices::Integer
    # all below have length nVertices, index i in any list is about vertex #i
    latVertex::AbstractArray{F} where F <: AbstractFloat
    lonVertex::AbstractArray{F} where F <: AbstractFloat
    xVertex::AbstractArray{F} where F <: AbstractFloat
    yVertex::AbstractArray{F} where F <: AbstractFloat
    vertexDegree::Integer
    cellsOnVertex::AbstractArray{I} where I <: Integer
    edgesOnVertex::AbstractArray{I} where I <: Integer
    edgeSignOnVertex::AbstractArray{I} where I <: Integer
    fVertex::AbstractArray{F} where F <: AbstractFloat
    areaTriangle::AbstractArray{F} where F <: AbstractFloat
    kiteAreasOnVertex::AbstractArray{F} where F <: AbstractFloat
    maxLevelVertexTop::AbstractArray{I} where I <: Integer
    maxLevelVertexBot::AbstractArray{I} where I <: Integer
    boundaryVertex::AbstractArray{I} where I <: Integer


    nVertLevels::Integer


    ExactSolutionParameters::AbstractArray{F} where F <: AbstractFloat
    myNamelist::Namelist

    gridSpacingMagnitude::AbstractFloat

    # X and Y size of grid, domain
    lX::AbstractFloat
    lY::AbstractFloat

    nNonPeriodicBoundaryCells::Integer
    nNonPeriodicBoundaryEdges::Integer
    nNonPeriodicBoundaryVertices::Integer



    device
    workGroupSize::Integer




    # load state from file
    function MPAS_Ocean(print_basic_geometry_mesh=false, mesh_directory = "MPAS_O_Shallow_Water/Mesh+Initial_Condition+Registry_Files/Periodic",
                        base_mesh_file_name = "base_mesh.nc",
                        mesh_file_name = "mesh.nc";
                        periodicity = "Periodic")
        myMPAS_O = new()




        myMPAS_O.myNamelist = Namelist()

        base_mesh_file = NCDatasets.Dataset("$mesh_directory/$base_mesh_file_name", "r", format=:netcdf4_classic)
        mesh_file = NCDatasets.Dataset("$mesh_directory/$mesh_file_name", "r", format=:netcdf4_classic)

        maxEdges = mesh_file.dim["maxEdges"]

        ## defining the mesh
        myMPAS_O.nCells = mesh_file.dim["nCells"]
        myMPAS_O.nEdges = mesh_file.dim["nEdges"]
        myMPAS_O.nVertices = mesh_file.dim["nVertices"]
        myMPAS_O.vertexDegree = mesh_file.dim["vertexDegree"]
        myMPAS_O.nVertLevels = 1

        # will be incremented in ocn_init_routines_compute_max_level
        myMPAS_O.nNonPeriodicBoundaryEdges = 0
        myMPAS_O.nNonPeriodicBoundaryVertices = 0
        myMPAS_O.nNonPeriodicBoundaryCells = 0

        nCells = myMPAS_O.nCells
        nEdges = myMPAS_O.nEdges
        nVertices = myMPAS_O.nVertices

        # Choose my_mesh_file_name to be either base_mesh_file_name or mesh_file_name
        my_mesh_file_name = mesh_file_name
        # Choose my_mesh_file to be either base_mesh_file or mesh_file
        my_mesh_file = mesh_file

        myMPAS_O.cellsOnCell = my_mesh_file["cellsOnCell"][:]
        myMPAS_O.nEdgesOnCell = base_mesh_file["nEdgesOnCell"][:]
        myMPAS_O.edgesOnCell = my_mesh_file["edgesOnCell"][:]
        myMPAS_O.verticesOnCell = my_mesh_file["verticesOnCell"][:]

        myMPAS_O.cellsOnEdge = my_mesh_file["cellsOnEdge"][:]
        myMPAS_O.nEdgesOnEdge = base_mesh_file["nEdgesOnEdge"][:]
        myMPAS_O.edgesOnEdge = my_mesh_file["edgesOnEdge"][:]
        myMPAS_O.verticesOnEdge = my_mesh_file["verticesOnEdge"][:]

        myMPAS_O.cellsOnVertex = my_mesh_file["cellsOnVertex"][:]
        myMPAS_O.edgesOnVertex = my_mesh_file["edgesOnVertex"][:]
        #


        myMPAS_O.angleEdge = my_mesh_file["angleEdge"][:]
        # if true
        myMPAS_O.angleEdge[:] = fix_angleEdge(mesh_directory,my_mesh_file_name,determineYCellAlongLatitude=true,
                                   printOutput=false,printRelevantMeshData=false)
        # end



        myMPAS_O.xEdge = my_mesh_file["xEdge"][:]
        myMPAS_O.yEdge = my_mesh_file["yEdge"][:]


        myMPAS_O.latCell = base_mesh_file["latCell"][:]
        myMPAS_O.lonCell = base_mesh_file["lonCell"][:]
        myMPAS_O.xCell = base_mesh_file["xCell"][:]
        myMPAS_O.yCell = base_mesh_file["yCell"][:]

        myMPAS_O.areaCell = my_mesh_file["areaCell"][:]
        myMPAS_O.dvEdge = my_mesh_file["dvEdge"][:]
        myMPAS_O.dcEdge = my_mesh_file["dcEdge"][:]

        myMPAS_O.kiteAreasOnVertex = my_mesh_file["kiteAreasOnVertex"]
        myMPAS_O.areaTriangle = my_mesh_file["areaTriangle"][:]


        myMPAS_O.latVertex = base_mesh_file["latVertex"][:]
        myMPAS_O.lonVertex = base_mesh_file["lonVertex"][:]
        myMPAS_O.xVertex = base_mesh_file["xVertex"][:]
        myMPAS_O.yVertex = base_mesh_file["yVertex"][:]


        myMPAS_O.gridSpacing = mesh_file["gridSpacing"][:]

        useGridSpacingMagnitudeDefaultDefinition = true
        if useGridSpacingMagnitudeDefaultDefinition
            myMPAS_O.gridSpacingMagnitude = myMPAS_O.gridSpacing[0+1]
        else
            myMPAS_O.gridSpacingMagnitude = myMPAS_O.xCell[1+1] - myMPAS_O.xCell[0+1]
        end


        # For a mesh with non-periodic zonal boundaries, adjust the zonal coordinates of the cell centers,
        # vertices and edges.
        dx = myMPAS_O.gridSpacingMagnitude # i.e. dx = myMPAS_O.dcEdge[0]
        dy = sqrt(3.0)/2.0*dx
        if periodicity == "NonPeriodic_x" || periodicity == "NonPeriodic_xy"
            myMPAS_O.xCell[:] .-= dx
            myMPAS_O.xVertex[:] .-= dx
            myMPAS_O.xEdge[:] .-= dx
        end
        if periodicity == "NonPeriodic_y" || periodicity == "NonPeriodic_xy"
            myMPAS_O.yCell[:] .-= dy
            myMPAS_O.yVertex[:] .-= dy
            myMPAS_O.yEdge[:] .-= dy
        end

        # Specify the zonal && meridional extents of the domain.
        myMPAS_O.lX = round(maximum(myMPAS_O.xCell)) # i.e. myMPAS_O.lX = sqrt(float(myMPAS_O.nCells))*dx
        # Please note that for all of our test problems, the MPAS-O mesh is generated in such a way that
        # myMPAS_O.lX i.e. the number of cells in the x (|| y) direction times dx has zero fractional part in
        # units of m, for which we can afford to round it to attain perfection.
        myMPAS_O.lY = sqrt(3.0)/2.0*myMPAS_O.lX # i.e. myMPAS_O.lY = max(myMPAS_O.yVertex)


        myMPAS_O.ExactSolutionParameters = DetermineExactSolutionParameters(myMPAS_O,false)
        # Define && initialize the following arrays not contained within either the base mesh file || the mesh
        # file.
        myMPAS_O.fVertex = zeros(Float64, nVertices)
        myMPAS_O.fCell = zeros(Float64, nCells)
        myMPAS_O.fEdge = zeros(Float64, nEdges)
        myMPAS_O.bottomDepth = zeros(Float64, nCells)
        DetermineCoriolisParameterAndBottomDepth(myMPAS_O)

        myMPAS_O.kiteIndexOnCell = zeros(Int64, (nCells,maxEdges))
        myMPAS_O.edgeSignOnVertex = zeros(Int8, (nVertices,maxEdges))


        myMPAS_O.weightsOnEdge = my_mesh_file["weightsOnEdge"][:]

        myMPAS_O.edgeSignOnCell = zeros(Int8, (nCells,maxEdges))
        ocn_init_routines_setup_sign_and_index_fields!(myMPAS_O)


        myMPAS_O.boundaryCell = zeros(Int8, (nCells, myMPAS_O.nVertLevels))
        myMPAS_O.boundaryEdge = zeros(Int8, (nEdges, myMPAS_O.nVertLevels))
        myMPAS_O.boundaryVertex = zeros(Int8, (nVertices, myMPAS_O.nVertLevels))
        myMPAS_O.edgeMask = zeros(Int8, nEdges)

        myMPAS_O.maxLevelCell = zeros(Int64, nCells)
        myMPAS_O.maxLevelEdgeTop = zeros(Int64, nEdges)
        myMPAS_O.maxLevelEdgeBot = zeros(Int64, nEdges)
        myMPAS_O.maxLevelVertexTop = zeros(Int64, nVertices)
        myMPAS_O.maxLevelVertexBot = zeros(Int64, nVertices)


        # replace!(myMPAS_O.edgesOnEdge, 0=>myMPAS_O.nEdges)
        # replace!(myMPAS_O.cellsOnEdge, 0=>myMPAS_O.nCells)


        ## defining the prognostic variables
        myMPAS_O.sshCurrent = zeros(Float64, nCells)
        myMPAS_O.sshNew = zeros(Float64, nCells)

        myMPAS_O.normalVelocityCurrent = zeros(Float64, (nEdges))
        myMPAS_O.normalVelocityNew = zeros(Float64, (nEdges))


        myMPAS_O.workGroupSize = 16

        return myMPAS_O
    end
end






function moveToDevice!(myMPAS_O, device)
    myMPAS_O.device = device

    arraytype = myMPAS_O.device == CUDAKernels.CUDADevice() ? CUDA.CuArray : Array

    for field in fieldnames(MPAS_Ocean)
        if typeof(getfield(myMPAS_O, field)) <: AbstractArray
            setfield!(myMPAS_O, field, arraytype(getfield(myMPAS_O, field)))
        end
    end
end





function ocn_init_routines_setup_sign_and_index_fields!(myMPAS_O)
    for iCell in range(1,myMPAS_O.nCells,step=1)
        for i in range(1,myMPAS_O.nEdgesOnCell[iCell],step=1)
            iEdge = myMPAS_O.edgesOnCell[i,iCell]
            iVertex = myMPAS_O.verticesOnCell[i,iCell]
            # Vector points from cell 1 to cell 2
            if myMPAS_O.cellsOnEdge[1,iEdge] == iCell
                myMPAS_O.edgeSignOnCell[iCell,i] = -1
            else
                myMPAS_O.edgeSignOnCell[iCell,i] = 1
            end
            for j in range(1,myMPAS_O.vertexDegree,step=1)
                if myMPAS_O.cellsOnVertex[j,iVertex] == iCell
                    myMPAS_O.kiteIndexOnCell[iCell,i] = j
                end
            end
        end
    end
#     printKiteIndexOnCell = false
#     if printKiteIndexOnCell
#         for iCell in range(1,myMPAS_O.nCells,step=1)
#             for i in range(1,myMPAS_O.nEdgesOnCell[iCell],step=1)
#                 println("On edge $i of cell $iCell, kiteIndexOnCell is $(myMPAS_O.kiteIndexOnCell[iCell,i]).")
#             end
#         end
#     end

    for iVertex in range(1,myMPAS_O.nVertices,step=1)
        for i in range(1,myMPAS_O.vertexDegree,step=1)
            iEdge = myMPAS_O.edgesOnVertex[i,iVertex]
            if iEdge !== 0
                # Vector points from vertex 1 to vertex 2
                if myMPAS_O.verticesOnEdge[1,iEdge] == iVertex
                    myMPAS_O.edgeSignOnVertex[iVertex,i] = -1
                else
                    myMPAS_O.edgeSignOnVertex[iVertex,i] = 1
                end
            end
        end
    end
end



function ocn_init_routines_compute_max_level!(myMPAS_O)
    for iEdge in 1:myMPAS_O.nEdges
        iCell1 = myMPAS_O.cellsOnEdge[1,iEdge]
        iCell2 = myMPAS_O.cellsOnEdge[2,iEdge]
        if iCell1 == 0 || iCell2 == 0
            myMPAS_O.boundaryEdge[iEdge,:] .= 1.0
            myMPAS_O.maxLevelEdgeTop[iEdge] = -1
            if iCell1 == 0
                myMPAS_O.maxLevelEdgeBot[iEdge] = myMPAS_O.maxLevelCell[iCell2]
            elseif iCell2 == 0
                myMPAS_O.maxLevelEdgeBot[iEdge] = myMPAS_O.maxLevelCell[iCell1]
            end
        elseif ! ( iCell1 == 0 || iCell2 == 0 )
            myMPAS_O.maxLevelEdgeTop[iEdge] = min(myMPAS_O.maxLevelCell[iCell1], myMPAS_O.maxLevelCell[iCell2])
            myMPAS_O.maxLevelEdgeBot[iEdge] = max(myMPAS_O.maxLevelCell[iCell1], myMPAS_O.maxLevelCell[iCell2])
        end
    end

    for iVertex in 1:myMPAS_O.nVertices
        iCell1 = myMPAS_O.cellsOnVertex[1,iVertex]
        if iCell1 == 0
            myMPAS_O.maxLevelVertexBot[iVertex] = -1
            myMPAS_O.maxLevelVertexTop[iVertex] = -1
        else
            myMPAS_O.maxLevelVertexBot[iVertex] = myMPAS_O.maxLevelCell[iCell1]
            myMPAS_O.maxLevelVertexTop[iVertex] = myMPAS_O.maxLevelCell[iCell1]
        end

        for i in 1:myMPAS_O.vertexDegree
            iCell = myMPAS_O.cellsOnVertex[i, iVertex]
            if iCell == 0
                myMPAS_O.maxLevelVertexBot[iVertex] = max(myMPAS_O.maxLevelVertexBot[iVertex], -1)
                myMPAS_O.maxLevelVertexBot[iVertex] = min(myMPAS_O.maxLevelVertexTop[iVertex], -1)
            else
                myMPAS_O.maxLevelVertexBot[iVertex] = max(myMPAS_O.maxLevelVertexBot[iVertex],
                                                            myMPAS_O.maxLevelCell[iCell])
                myMPAS_O.maxLevelVertexTop[iVertex] = min(myMPAS_O.maxLevelVertexTop[iVertex],
                                                            myMPAS_O.maxLevelCell[iCell])
            end
        end
    end

    determine_boundaryEdge_Generalized_Method = true

    if determine_boundaryEdge_Generalized_Method
        myMPAS_O.boundaryEdge[:,:] .= 1
    end
    myMPAS_O.edgeMask[:,:] .= 0

    for iEdge in 1:myMPAS_O.nEdges
        index_UpperLimit = myMPAS_O.maxLevelEdgeTop[iEdge]
        if index_UpperLimit > -1
            if determine_boundaryEdge_Generalized_Method
                myMPAS_O.boundaryEdge[iEdge,1:index_UpperLimit+1] .= 0
            end
            myMPAS_O.edgeMask[iEdge,1:index_UpperLimit+1] .= 1
        end
    end

    for iEdge in 1:myMPAS_O.nEdges
        iCell1 = myMPAS_O.cellsOnEdge[1, iEdge]
        iCell2 = myMPAS_O.cellsOnEdge[2, iEdge]

        iVertex1 = myMPAS_O.verticesOnEdge[1, iEdge]
        iVertex2 = myMPAS_O.verticesOnEdge[2, iEdge]

        if myMPAS_O.boundaryEdge[iEdge] == 1
            if iCell1 !== 0
                myMPAS_O.boundaryCell[iCell1] = 1
            end
            if iCell2 !== 0
                myMPAS_O.boundaryCell[iCell2] = 1
            end
            myMPAS_O.boundaryVertex[iVertex1] = 1
            myMPAS_O.boundaryVertex[iVertex2] = 1
        end
    end


    for iCell in 1:myMPAS_O.nCells
        for k in 1:myMPAS_O.nVertLevels
            if myMPAS_O.maxLevelCell[iCell] >= k
                myMPAS_O.cellMask[iCell, k] = 1
            end
        end
    end
    for iVertex in 1:myMPAS_O.nVertices
        for k in 1:myMPAS_O.nVertLevels
            if myMPAS_O.maxLevelVertexBot[iVertex] >= k
                myMPAS_O.vertexMask[iVertex, k] = 1
            end
        end
    end


    for iEdge in 1:myMPAS_O.nEdges
        if myMPAS_O.boundaryEdge[iEdge, 1] == 1.0
            myMPAS_O.nNonPeriodicBoundaryEdges += 1
        end
    end
    for iVertex in 1:myMPAS_O.nVertices
        if myMPAS_O.boundaryVertex[iVertex, 1] == 1.0
            myMPAS_O.nNonPeriodicBoundaryVertices += 1
        end
    end
    for iCell in 1:myMPAS_O.nCells
        if myMPAS_O.boundaryCell[iCell, 1] == 1.0
            myMPAS_O.nNonPeriodicBoundaryCells += 1
        end
    end
end


function DetermineExactSolutionParameters(myMPAS_O,printPhaseSpeedOfWaveModes)
    angleEdge = 0.0
    beta0 = myMPAS_O.myNamelist.config_meridional_gradient_of_Coriolis_parameter
    c = myMPAS_O.myNamelist.config_phase_speed_of_coastal_Kelvin_wave
    dx = myMPAS_O.gridSpacingMagnitude # i.e. dx = myMPAS_O.dcEdge[0+1]
    dy = sqrt(3.0)/2.0*dx
    etaHat1 = myMPAS_O.myNamelist.config_surface_elevation_amplitude
    if myMPAS_O.myNamelist.config_problem_type == "Viscous_Burgers_Equation"
        etaHat2 = 0.0
    else
        etaHat2 = 2.0*etaHat1
    end
    f0 = myMPAS_O.myNamelist.config_Coriolis_parameter
    g = myMPAS_O.myNamelist.config_gravity
    H = myMPAS_O.myNamelist.config_mean_depth
    if myMPAS_O.myNamelist.config_problem_type == "Topographic_Rossby_Wave"
        alpha0 = beta0*H/f0
        # In the Northern Hemisphere where f0 > 0, the topographic Rossby wave travels with the shallower water on
        # its right. Hence if alpha0 > 0 i.e. the ocean depth increases northward, the topographic Rossby wave
        # will propagate eastward else it will propagate westward.
        myMPAS_O.myNamelist.config_bottom_slope = alpha0
    else
        alpha0 = myMPAS_O.myNamelist.config_bottom_slope # alpha0 = 0.0
    end
    if myMPAS_O.myNamelist.config_problem_type == "Barotropic_Tide"
        kX1 = 2.5*pi/myMPAS_O.lX
        kY1 = 0.0
        kX2 = 4.5*pi/myMPAS_O.lX
        kY2 = 0.0
    elseif myMPAS_O.myNamelist.config_problem_type == "Viscous_Burgers_Equation"
        kX1 = 0.0
        kY1 = 0.0
        kX2 = 0.0
        kY2 = 0.0
    else
        kX1 = 2.0*pi/myMPAS_O.lX
        kY1 = 2.0*pi/myMPAS_O.lY
        kX2 = 2.0*kX1
        kY2 = 2.0*kY1
    end
    lX = myMPAS_O.lX
    lY = myMPAS_O.lY
    R = myMPAS_O.myNamelist.config_radius_of_deformation
    Req = myMPAS_O.myNamelist.config_equatorial_radius_of_deformation
    LengthScale = sqrt(c/beta0)
    TimeScale = 1.0/sqrt(beta0*c)
    VelocityScale = c
    SurfaceElevationScale = c^2.0/g
    if (myMPAS_O.myNamelist.config_problem_type == "default"
        || myMPAS_O.myNamelist.config_problem_type == "Coastal_Kelvin_Wave")
        omega1 = -c*kY1
        omega2 = -c*kY2
    elseif myMPAS_O.myNamelist.config_problem_type == "Inertia_Gravity_Wave"
        omega1 = sqrt(g*H*(kX1^2.0 + kY1^2.0) + f0^2.0)
        omega2 = sqrt(g*H*(kX2^2.0 + kY2^2.0) + f0^2.0)
    elseif myMPAS_O.myNamelist.config_problem_type == "Planetary_Rossby_Wave"
        omega1 = -beta0*R^2.0*kX1/(1.0 + R^2.0*(kX1^2.0 + kY1^2.0))
        omega2 = -beta0*R^2.0*kX2/(1.0 + R^2.0*(kX2^2.0 + kY2^2.0))
    elseif myMPAS_O.myNamelist.config_problem_type == "Topographic_Rossby_Wave"
        omega1 = alpha0*g/f0*kX1/(1.0 + R^2.0*(kX1^2.0 + kY1^2.0))
        omega2 = alpha0*g/f0*kX2/(1.0 + R^2.0*(kX2^2.0 + kY2^2.0))
    elseif myMPAS_O.myNamelist.config_problem_type == "Equatorial_Kelvin_Wave"
        omega1 = c*kX1
        omega2 = c*kX2
    elseif myMPAS_O.myNamelist.config_problem_type == "Equatorial_Yanai_Wave"
        omega1 = GWESST.DetermineEquatorialYanaiWaveNonDimensionalAngularFrequency(LengthScale*kX1)/TimeScale
        omega2 = GWESST.DetermineEquatorialYanaiWaveNonDimensionalAngularFrequency(LengthScale*kX2)/TimeScale
    elseif myMPAS_O.myNamelist.config_problem_type == "Equatorial_Rossby_Wave"
        omega1 = GWESST.DetermineEquatorialRossbyAndInertiaGravityWaveNonDimensionalAngularFrequency(
                 "Equatorial_Rossby_Wave",LengthScale*kX1,m=1)/TimeScale
        omega2 = GWESST.DetermineEquatorialRossbyAndInertiaGravityWaveNonDimensionalAngularFrequency(
                 "Equatorial_Rossby_Wave",LengthScale*kX2,m=1)/TimeScale
    elseif myMPAS_O.myNamelist.config_problem_type == "Equatorial_Inertia_Gravity_Wave"
        omega1 = GWESST.DetermineEquatorialRossbyAndInertiaGravityWaveNonDimensionalAngularFrequency(
                 "Equatorial_Inertia_Gravity_Wave",LengthScale*kX1,m=2)/TimeScale
        omega2 = GWESST.DetermineEquatorialRossbyAndInertiaGravityWaveNonDimensionalAngularFrequency(
                 "Equatorial_Inertia_Gravity_Wave",LengthScale*kX2,m=2)/TimeScale
    elseif myMPAS_O.myNamelist.config_problem_type == "Barotropic_Tide"
        omega1 = sqrt(g*H*kX1^2.0 + f0^2.0)
        omega2 = sqrt(g*H*kX2^2.0 + f0^2.0)
    elseif (myMPAS_O.myNamelist.config_problem_type == "Diffusion_Equation"
          || myMPAS_O.myNamelist.config_problem_type == "Viscous_Burgers_Equation")
        omega1 = 0.0
        omega2 = 0.0
    end
    if (myMPAS_O.myNamelist.config_problem_type == "default"
        || myMPAS_O.myNamelist.config_problem_type == "Coastal_Kelvin_Wave")
        cX1 = 0.0
        cY1 = -c # cY1 = omega1/kY1 = -c
        cX2 = 0.0
        cY2 = -c # cY2 = omega2/kY2 = -c
    elseif myMPAS_O.myNamelist.config_problem_type_Equatorial_Wave
        if myMPAS_O.myNamelist.config_problem_type == "Equatorial_Kelvin_Wave"
            cX1 = c
            cY1 = 0.0
            cX2 = c
            cY2 = 0.0
        else
            cX1 = omega1/kX1
            cY1 = 0.0
            cX2 = omega2/kX2
            cY2 = 0.0
        end
    elseif myMPAS_O.myNamelist.config_problem_type == "Barotropic_Tide"
        cX1 = omega1/kX1
        cY1 = 0.0
        cX2 = omega2/kX2
        cY2 = 0.0
    elseif (myMPAS_O.myNamelist.config_problem_type == "Diffusion_Equation"
          || myMPAS_O.myNamelist.config_problem_type == "Viscous_Burgers_Equation")
        cX1 = 0.0
        cY1 = 0.0
        cX2 = 0.0
        cY2 = 0.0
    else
        cX1 = omega1/kX1
        cY1 = omega1/kY1
        cX2 = omega2/kX2
        cY2 = omega2/kY2
    end
    if ((myMPAS_O.myNamelist.config_problem_type_Geophysical_Wave
         || myMPAS_O.myNamelist.config_problem_type == "Barotropic_Tide") && printPhaseSpeedOfWaveModes)
        println("The zonal component of the phase speed of the first wave mode is %.4g." %cX1)
        println("The meridional component of the phase speed of the first wave mode is %.4g." %cY1)
        println("The zonal component of the phase speed of the second wave mode is %.4g." %cX2)
        println("The meridional component of the phase speed of the second wave mode is %.4g." %cY2)
    end
    kappaX = myMPAS_O.myNamelist.config_zonal_diffusivity
    kappaY = myMPAS_O.myNamelist.config_meridional_diffusivity
    kappa1 = kappaX*kX1^2.0 + kappaY*kY1^2.0
    kappa2 = kappaX*kX2^2.0 + kappaY*kY2^2.0
    uL = myMPAS_O.myNamelist.config_viscous_Burgers_zonal_velocity_left
    uR = myMPAS_O.myNamelist.config_viscous_Burgers_zonal_velocity_right
    s = myMPAS_O.myNamelist.config_viscous_Burgers_shock_speed
    myMPAS_O.myNamelist.config_viscous_Burgers_zonal_offset = 0.25*lX
    x0 = myMPAS_O.myNamelist.config_viscous_Burgers_zonal_offset
    ExactSolutionParameters = [alpha0,angleEdge,beta0,c,cX1,cX2,cY1,cY2,dx,dy,etaHat1,etaHat2,f0,g,H,kX1,kX2,kY1,
                               kY2,lX,lY,omega1,omega2,R,Req,LengthScale,TimeScale,VelocityScale,
                               SurfaceElevationScale,kappaX,kappaY,kappa1,kappa2,uL,uR,s,x0]
    return ExactSolutionParameters
end



function DetermineCoriolisParameterAndBottomDepth(myMPAS_O)
    alpha0 = myMPAS_O.ExactSolutionParameters[0+1]
    beta0 = myMPAS_O.ExactSolutionParameters[2+1]
    f0 = myMPAS_O.ExactSolutionParameters[12+1]
    H = myMPAS_O.ExactSolutionParameters[14+1]
    if (myMPAS_O.myNamelist.config_problem_type == "default"
        || myMPAS_O.myNamelist.config_problem_type == "Barotropic_Tide"
        || myMPAS_O.myNamelist.config_problem_type_Geophysical_Wave)
        if myMPAS_O.myNamelist.config_problem_type == "Planetary_Rossby_Wave"
            myMPAS_O.fCell[:] = f0 + beta0*myMPAS_O.yCell[:]
            myMPAS_O.fEdge[:] = f0 + beta0*myMPAS_O.yEdge[:]
            myMPAS_O.fVertex[:] = f0 + beta0*myMPAS_O.yVertex[:]
        else
            myMPAS_O.fCell[:] .= f0
            myMPAS_O.fEdge[:] .= f0
            myMPAS_O.fVertex[:] .= f0
        end
        if myMPAS_O.myNamelist.config_problem_type == "Topographic_Rossby_Wave"
            myMPAS_O.bottomDepth[:] = H + alpha0*yCell[:]
        else
            myMPAS_O.bottomDepth[:] .= H
        end
    end
end
