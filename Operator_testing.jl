#!/usr/bin/env python
# coding: utf-8

# # Spatial Operator Testing
# #### Numerical methods are compared to exact analytical solutions

# # Table of Contents <a id="toc"></a>
# Convergence Tests:
#  * [Tangential Velocity](#tangentialvelocity)
#  * Curl
#      * [vertex](#curlvertex)
#      * [cell center](#curlcellcenter)
#  * [Gradient](#gradient)
#  * [Divergence](#divergence)

# In[1]:


CODE_ROOT = pwd() * "/"


# In[2]:


include(CODE_ROOT * "mode_forward/time_steppers.jl")
include(CODE_ROOT * "mode_init/MPAS_Ocean.jl")
include(CODE_ROOT * "visualization.jl")


using PyPlot
using PyCall

animation  = pyimport("matplotlib.animation")
ipydisplay = pyimport("IPython.display")


using LinearAlgebra # for norm()

using Printf # for print formatting


# In[3]:


using DelimitedFiles
using Dates


# In[4]:


gravity = 9.8


# In[5]:


device = "CPU"


# ## define continuous ocean and analytical form of operators in the form of functions
# functions such as surface_elevation_1 calculate and return the value of the surface elevation at x and y

# In[6]:


function surface_elevation_1(lX, lY, x, y)
    eta0 = 0.1
    eta = eta0 * sin( 2.0 * pi * x / lX ) * sin( 2.0 * pi * y / lY )
    return eta
end


# In[7]:


function problem_specific_prefix_1()
    prefix = "Expt1_"
    return prefix
end


# In[8]:


function surface_elevation_gradient_1(lX, lY, x, y)
    eta0 = 0.1
    eta = eta0 * sin( 2.0 * pi * x / lX ) * sin( 2.0 * pi * y / lX )
    
    eta_x = eta0 * (2.0 * pi / lX)*cos( 2.0 * pi * x /lX ) * sin( 2.0 * pi * y/lY)
    eta_y = eta0 * (2.0 * pi / lY)*cos( 2.0 * pi * y /lY ) * sin( 2.0 * pi * x/lX)
    return eta_x, eta_y
end


# In[9]:


function surface_elevation_laplacian_1(lX,lY,x,y)
    eta0 = 0.1
    eta = eta0 * sin( 2.0 * pi * x / lX ) * sin( 2.0 * pi * y / lY )
    eta_xx = - ( 2.0 * pi / lX )^2.0 * eta
    eta_yy = - ( 2.0 * pi / lY )^2.0 * eta
    laplacian = eta_xx + eta_yy
    return laplacian
end


# In[10]:


surface_elevation = surface_elevation_1

surface_elevation_gradient = surface_elevation_gradient_1

surface_elevation_laplacian = surface_elevation_laplacian_1

problem_specific_prefix = problem_specific_prefix_1


# In[11]:


function velocity(lX,lY,x,y, f=1e-4, g=10.0)
    eta_x, eta_y = surface_elevation_gradient(lX,lY,x,y) 
    v = -gravity*eta_x/f
    u = gravity*eta_y/f
    return u, v
end


# In[12]:


function velocity_curl(lX,lY,x,y, f=1e-4, g=10.0)
    zeta = gravity/f*surface_elevation_laplacian(lX,lY,x,y)
    return zeta
end


# In[13]:


function ComputeNormalAndTangentialComponentsAtEdge(myVectorQuantityAtEdge,angleEdge,returningComponent)
    nEdges = length(angleEdge)
    if returningComponent == "normal" || returningComponent == "both"
        myVectorQuantityAtEdgeNormalComponent = zeros(nEdges)
    end
    if returningComponent == "tangential" || returningComponent == "both"
        myVectorQuantityAtEdgeTangentialComponent = zeros(nEdges)
    end
    for iEdge in range(0+1,nEdges,step=1)
        xComponent = myVectorQuantityAtEdge[iEdge,0+1]
        yComponent = myVectorQuantityAtEdge[iEdge,1+1]
        if returningComponent == "normal" || returningComponent == "both"
            myVectorQuantityAtEdgeNormalComponent[iEdge] = (xComponent*cos(angleEdge[iEdge]) 
                                                            + yComponent*sin(angleEdge[iEdge]))
        end
        if returningComponent == "tangential" || returningComponent == "both"
            myVectorQuantityAtEdgeTangentialComponent[iEdge] = (yComponent*cos(angleEdge[iEdge]) 
                                                            - xComponent*sin(angleEdge[iEdge]))
        end
    end
    if returningComponent == "normal"
        return myVectorQuantityAtEdgeNormalComponent
    elseif returningComponent == "tangential"
        return myVectorQuantityAtEdgeTangentialComponent
    else # if returningComponent == "both"
        return myVectorQuantityAtEdgeNormalComponent, myVectorQuantityAtEdgeTangentialComponent
    end
end


# # Tangential Velocity Test

# ## define numerical method for tangential velocity

# In[14]:


function numerical_tangential_velocity(mpasOcean,myNormalVelocity,periodicity)
    myTangentialVelocity = zeros(mpasOcean.nEdges)
    for iEdge in range(0+1,mpasOcean.nEdges,step=1)
        if ((periodicity == "NonPeriodic_x" || periodicity == "NonPeriodic_y" || periodicity == "NonPeriodic_xy") 
            && mpasOcean.boundaryEdge[iEdge] == 1.0)
            u, v = velocity(mpasOcean.lX,mpasOcean.lY,mpasOcean.xEdge[iEdge],mpasOcean.yEdge[iEdge])
            myTangentialVelocity[iEdge] = v*cos(mpasOcean.angleEdge[iEdge]) - u*sin(mpasOcean.angleEdge[iEdge])
        else
            myTangentialVelocity[iEdge] = 0.0
            # Compute tangential velocities
            for i in range(0+1,mpasOcean.nEdgesOnEdge[iEdge],step=1)
                eoe = mpasOcean.edgesOnEdge[i,iEdge]
                weightOnEdge = mpasOcean.weightsOnEdge[i,iEdge]
                myTangentialVelocity[iEdge] += weightOnEdge*myNormalVelocity[eoe]
            end
        end
    end
    return myTangentialVelocity
end


# ## compare numerical and analytical tangential velocity

# In[15]:


function test_tangential_velocity(plotFigures,mesh_directory,base_mesh_file_name,mesh_file_name,periodicity)
    
    mpasOcean = MPAS_Ocean(mesh_directory,base_mesh_file_name,mesh_file_name,
                               periodicity=periodicity)
    
#     moveArrays!(mpasOcean, arraytype)
        
    if periodicity == "NonPeriodic_x"
        iEdgeStartingIndex = 1+1
    else
        iEdgeStartingIndex = 0+1
    end
    
    prefix = problem_specific_prefix()
    myAnalyticalVelocityComponentsAtEdge = zeros((mpasOcean.nEdges,2))
    myAnalyticalNormalVelocity = zeros(mpasOcean.nEdges)
    myAnalyticalTangentialVelocity = zeros(mpasOcean.nEdges)
    
    for iEdge in range(0+1,mpasOcean.nEdges,step=1)
        myAnalyticalVelocityComponentsAtEdge[iEdge,0+1], myAnalyticalVelocityComponentsAtEdge[iEdge,1+1] = (
            velocity(mpasOcean.lX,mpasOcean.lY,mpasOcean.xEdge[iEdge],mpasOcean.yEdge[iEdge])
        )
    end
    
#     fig, _ = edgeHeatMapMesh(mpasOcean, myAnalyticalVelocityComponentsAtEdge[:,1])
#     display(fig)
#     fig, _ = edgeHeatMapMesh(mpasOcean, myAnalyticalVelocityComponentsAtEdge[:,2])
#     display(fig)
    
    myAnalyticalNormalVelocity, myAnalyticalTangentialVelocity = (
                ComputeNormalAndTangentialComponentsAtEdge(myAnalyticalVelocityComponentsAtEdge,mpasOcean.angleEdge,"both"))
    myNumericalTangentialVelocity = (
                numerical_tangential_velocity(mpasOcean,myAnalyticalNormalVelocity,periodicity))
    myTangentialVelocityError = myNumericalTangentialVelocity - myAnalyticalTangentialVelocity
    MaxErrorNorm = norm(myTangentialVelocityError,Inf)
    L2ErrorNorm = (norm(myTangentialVelocityError)
                   /sqrt(float(mpasOcean.nEdges - mpasOcean.nNonPeriodicBoundaryEdges)))
    @sprintf("The maximum error norm of the tangential velocity is %.3g.", MaxErrorNorm)
    @sprintf("The L2 error norm of the tangential velocity is %.3g.", L2ErrorNorm)  
    
    if plotFigures
        xLabel = "Zonal Distance (km)"
        yLabel = "Meridional Distance (km)"
        
        
        Title = "Analytical Tangential Velocity"
        FigureTitle = prefix * "TangentialVelocity_Analytical"
        fig, ax, _, _ = edgeHeatMapMesh(mpasOcean, myAnalyticalTangentialVelocity)
        ax.set_xlabel(xLabel)
        ax.set_ylabel(yLabel)
        ax.set_title(Title)
        fig.suptitle(FigureTitle)
        display(fig)
        
        Title = "Numerical Tangential Velocity"
        FigureTitle = prefix * "TangentialVelocity_Numerical"
        fig, ax, _, _  = edgeHeatMapMesh(mpasOcean, myNumericalTangentialVelocity)
        ax.set_xlabel(xLabel)
        ax.set_ylabel(yLabel)
        ax.set_title(Title)
        fig.suptitle(FigureTitle)
        display(fig)
        
        Title = "Tangential Velocity Error"
        FigureTitle = prefix * "TangentialVelocity_Error"
        fig, ax, _, _ = edgeHeatMapMesh(mpasOcean, myTangentialVelocityError)
        ax.set_xlabel(xLabel)
        ax.set_ylabel(yLabel)
        ax.set_title(Title)
        fig.suptitle(FigureTitle)
        display(fig)
    end
    
    return mpasOcean.nCells, MaxErrorNorm, L2ErrorNorm
end


# In[16]:


include(CODE_ROOT * "mode_init/MPAS_Ocean.jl")

# ## run convergence test to see how error scales with cell width

# In[18]:


function convergence_test(periodicity, mesh_directory, operator_name, test; resolutions=[64, 96, 144, 216, 324], showplots=true, savedata=false)
    nCases = length(resolutions)
    nCellsX = collect(Int.(round.(resolutions)))
    ncells = zeros(Float64, nCases)
    dts = zeros(Float64, nCases)
    MaxErrorNorm = zeros(Float64, nCases)
    L2ErrorNorm = zeros(Float64, nCases)
    
    for iCase = 1:nCases
        if periodicity == "Periodic"
            base_mesh_file_name = "base_mesh_$(nCellsX[iCase])x$(nCellsX[iCase]).nc"
        else
            base_mesh_file_name = "culled_mesh_$(nCellsX[iCase])x$(nCellsX[iCase]).nc"
        end
        mesh_file_name = "mesh_$(nCellsX[iCase])x$(nCellsX[iCase]).nc"
        println()
        println("running test $iCase of $nCases, mesh: $mesh_file_name")
        etaHat = 1e0
        ncells[iCase], MaxErrorNorm[iCase], L2ErrorNorm[iCase] =
                test(false, mesh_directory, base_mesh_file_name, mesh_file_name, periodicity)
    end
    
    
#     nCases = 6
#     ncells = zeros(Int64, nCases)
#     MaxErrorNorm = zeros(Float64, nCases)
#     L2ErrorNorm = zeros(Float64, nCases)
    
#     for iCase = 1:nCases
#         if periodicity == "Periodic"
#             base_mesh_file_name = "base_mesh_$iCase.nc"
#         else
#             base_mesh_file_name = "culled_mesh_$iCase.nc"
#         end
#         mesh_file_name = "mesh_$iCase.nc"
#         ncells[iCase], MaxErrorNorm[iCase], L2ErrorNorm[iCase] =
#                 test(false, mesh_directory, base_mesh_file_name, mesh_file_name, periodicity)
#     end
    
    A = [log10.(ncells)    ones(length(ncells))]
    m, c = A \ log10.(MaxErrorNorm)
    y = m*log10.(ncells) .+ c
    y = 10 .^ y
    
    ymax = y
    
    if showplots
        fig, ax = subplots(1,1)
        loglog(ncells, MaxErrorNorm, label="Maximum Error Norm", marker="s", linestyle="None", color="black")
        loglog(ncells, y, label="Best Fit Line, slope=$(round(m,digits=2))")
        ax.set_title("Convergence of $operator_name", fontweight="bold")
        ax.legend(loc="upper right")
        ax.set_xlabel("Number of cells")
        ax.set_ylabel("Maximum Error Norm of $operator_name")
        grid(which="both")
    end
    
    
    A = [log10.(ncells)    ones(length(ncells))]
    m, c = A \ log10.(L2ErrorNorm)
    y = m*log10.(ncells) .+ c
    y = 10 .^ y
    
    if savedata
        fpath = CODE_ROOT * "output/operator_convergence/$operator_name/$periodicity/"
        mkpath(fpath)
        fname = "$fpath$(Dates.now()).txt"
        open(fname, "w") do io
            writedlm(io, [ncells, y, ymax, L2ErrorNorm, MaxErrorNorm])
        end
        println("saved to $fname")
    end
    
    if showplots
        fig, ax = subplots(1,1)
        loglog(ncells, L2ErrorNorm, label="\$L^2\$ Error Norm", marker="s", linestyle="None", color="black")
        loglog(ncells, y, label="Best Fit Line, slope=$(round(m,digits=2))")
        ax.set_title("Convergence of $operator_name", fontweight="bold")
        ax.legend(loc="upper right")
        ax.set_xlabel("Number of cells")
        ax.set_ylabel("\$L^2\$ Error Norm of $operator_name")
        grid(which="both")
    end
end


# <a id="tangentialvelocity"></a>
# [Back to Top](#toc)


# ## success
# plots have a slope ~2, meaning the error scales with the cell width like it should!

# In[19]:


function numerical_curl_operator_at_vertex(mpasOcean, myNormalVelocity, periodicity)
    if periodicity == "NonPeriodic_x" || periodicity == "NonPeriodic_y" || periodicity == "NonPeriodic_xy"
        ocn_init_routines_compute_max_level!(mpasOcean)
    end
    
    myCurlAtVertex = zeros(Float64, mpasOcean.nVertices)
    for iVertex = 1:mpasOcean.nVertices
        if (periodicity !== "Periodic" && mpasOcean.boundaryVertex[iVertex] == 1.0)
            myCurlAtVertex[iVertex] = velocity_curl(mpasOcean.lX, mpasOcean.lY, mpasOcean.xVertex[iVertex], mpasOcean.yVertex[iVertex])
        else
            for i = 1:mpasOcean.vertexDegree
                iEdge = mpasOcean.edgesOnVertex[i, iVertex]
                if iEdge == 0
                    println(mpasOcean.boundaryVertex[iVertex])
                end
                myCurlAtVertex[iVertex] += mpasOcean.dcEdge[iEdge] * myNormalVelocity[iEdge] * mpasOcean.edgeSignOnVertex[iVertex,i]
            end
            myCurlAtVertex[iVertex] /= mpasOcean.areaTriangle[iVertex]
        end
    end
    return myCurlAtVertex
end


# In[20]:


function test_curl_at_vertex(plotFigures,
                        mesh_directory, base_mesh_file_name, mesh_file_name, periodicity)
    

    mpasOcean = MPAS_Ocean(mesh_directory,base_mesh_file_name,mesh_file_name,
                                       periodicity=periodicity, nvlevels=1)

    
    myAnalyticalVelocity = zeros(Float64, (mpasOcean.nEdges, 2))
    for iEdge = 1:mpasOcean.nEdges
        myAnalyticalVelocity[iEdge,1], myAnalyticalVelocity[iEdge,2] = velocity(mpasOcean.lX, mpasOcean.lY,
                                                                                    mpasOcean.xEdge[iEdge], mpasOcean.yEdge[iEdge])
    end
    
    myAnalyticalNormalVelocity = ComputeNormalAndTangentialComponentsAtEdge(myAnalyticalVelocity, mpasOcean.angleEdge, "normal")
    
    
    myNumericalCurl = numerical_curl_operator_at_vertex(mpasOcean, myAnalyticalNormalVelocity, periodicity)
    
    
    myAnalyticalCurl = zeros(Float64, mpasOcean.nVertices)
    for iVertex = 1:mpasOcean.nVertices
        myAnalyticalCurl[iVertex] = velocity_curl(mpasOcean.lX, mpasOcean.lY, mpasOcean.xVertex[iVertex], mpasOcean.yVertex[iVertex])
    end
    
    
    myCurlError = myNumericalCurl - myAnalyticalCurl
    MaxErrorNorm = norm(myCurlError, Inf)
    L2ErrorNorm = norm(myCurlError) / sqrt(float(mpasOcean.nVertices - mpasOcean.nNonPeriodicBoundaryVertices))
    
    if plotFigures
        fig, ax, _, _ = vertexHeatMapMesh(mpasOcean, myAnalyticalCurl, dotsize=10)
        ax.set_title("Analytical Curl")
        ax.set_xlabel("Zonal Distance (km)")
        ax.set_ylabel("Meridional Distance (km)")
        display(fig)
        
        fig, ax, _, _ = vertexHeatMapMesh(mpasOcean, myNumericalCurl, dotsize=10)
        ax.set_title("Numerical Curl")
        ax.set_xlabel("Zonal Distance (km)")
        ax.set_ylabel("Meridional Distance (km)")
        display(fig)
        
        fig, ax, _, _ = vertexHeatMapMesh(mpasOcean, myCurlError, dotsize=10)
        ax.set_title("Curl Error")
        ax.set_xlabel("Zonal Distance (km)")
        ax.set_ylabel("Meridional Distance (km)")
        display(fig)
    end
    
    return mpasOcean.nCells, MaxErrorNorm, L2ErrorNorm
end

# ## now do test of curl at cell centers instead of vertices

# In[22]:


function numerical_curl_operator_at_cell_center(mpasOcean, myNormalVelocity, periodicity)
    myCurlAtVertex = numerical_curl_operator_at_vertex(mpasOcean, myNormalVelocity, periodicity)
    myCurlAtCellCenter = zeros(Float64, mpasOcean.nCells)
    for iCell = 1:mpasOcean.nCells
        for i = 1:mpasOcean.nEdgesOnCell[iCell]
            iVertex = mpasOcean.verticesOnCell[i, iCell]
            iKite = mpasOcean.kiteIndexOnCell[iCell, i] #- 1
            myCurlAtCellCenter[iCell] += mpasOcean.kiteAreasOnVertex[iKite, iVertex] * myCurlAtVertex[iVertex]
        end
        myCurlAtCellCenter[iCell] /= mpasOcean.areaCell[iCell]
    end
    return myCurlAtCellCenter
end


# In[23]:


function test_curl_at_cell_center(plotFigures,
                            mesh_directory, base_mesh_file_name, mesh_file_name, periodicity)
    
    mpasOcean = MPAS_Ocean(mesh_directory,base_mesh_file_name,mesh_file_name,
                                       periodicity=periodicity, nvlevels=1)

    
    myAnalyticalVelocity = zeros(Float64, (mpasOcean.nEdges, 2))
    for iEdge = 1:mpasOcean.nEdges
        myAnalyticalVelocity[iEdge,1], myAnalyticalVelocity[iEdge,2] = velocity(mpasOcean.lX, mpasOcean.lY, mpasOcean.xEdge[iEdge], mpasOcean.yEdge[iEdge])
    end
    
    myAnalyticalNormalVelocity = ComputeNormalAndTangentialComponentsAtEdge(myAnalyticalVelocity, mpasOcean.angleEdge, "normal")
    
    
    myNumericalCurl = numerical_curl_operator_at_cell_center(mpasOcean, myAnalyticalNormalVelocity, periodicity)
    
    
    myAnalyticalCurl = zeros(Float64, mpasOcean.nCells)
    for iCell = 1:mpasOcean.nCells
        myAnalyticalCurl[iCell] = velocity_curl(mpasOcean.lX, mpasOcean.lY, mpasOcean.xCell[iCell], mpasOcean.yCell[iCell])
    end
    
    
    myCurlError = myNumericalCurl - myAnalyticalCurl
    MaxErrorNorm = norm(myCurlError, Inf)
    L2ErrorNorm = norm(myCurlError) / sqrt(float(mpasOcean.nVertices - mpasOcean.nNonPeriodicBoundaryVertices))
    
    if plotFigures
        fig, ax, _, _ = heatMapMesh(mpasOcean, myAnalyticalCurl)
        ax.set_title("Analytical Curl")
        ax.set_xlabel("Zonal Distance (km)")
        ax.set_ylabel("Meridional Distance (km)")
        display(fig)
        
        fig, ax, _, _ = heatMapMesh(mpasOcean, myNumericalCurl)
        ax.set_title("Numerical Curl")
        ax.set_xlabel("Zonal Distance (km)")
        ax.set_ylabel("Meridional Distance (km)")
        display(fig)
        
        fig, ax, _, _ = heatMapMesh(mpasOcean, myCurlError)
        ax.set_title("Curl Error")
        ax.set_xlabel("Zonal Distance (km)")
        ax.set_ylabel("Meridional Distance (km)")
        display(fig)
    end
    
    return mpasOcean.nCells, MaxErrorNorm, L2ErrorNorm
end



# <a id="curlvertex"></a>
# [Back to Top](#toc)

# ## interesting detail:
# the curl error at the vertex scales linearly, with a slope of 1 on the log plot.
# but the curl error calculated at the cell centers scales quadratically, with a slope of 2 on the log plot.

# <a id="curlcellcenter"></a>
# [Back to Top](#toc)

# # gradient operator test

# In[25]:


function numerical_gradient_operator(mpasOcean, myScalarAtCellCenters, periodicity)
    if periodicity == "NonPeriodic_x" || periodicity == "NonPeriodic_y" || periodicity == "NonPeriodic_xy"
        ocn_init_routines_compute_max_level!(mpasOcean)
    end
    
    myNumericalGradientNormalToEdge = zeros(Float64, mpasOcean.nEdges)
    for iEdge = 1:mpasOcean.nEdges
        if ! (periodicity !== "Periodic" && mpasOcean.boundaryEdge[iEdge] == 1.0)
            cell1, cell2 = mpasOcean.cellsOnEdge[:,iEdge]
            myNumericalGradientNormalToEdge[iEdge] = ( myScalarAtCellCenters[cell2] - myScalarAtCellCenters[cell1] ) / mpasOcean.dcEdge[iEdge]
        end
    end
    
    return myNumericalGradientNormalToEdge
end


# In[26]:


function test_numerical_gradient_operator(plotFigures,
                                mesh_directory, base_mesh_file_name, mesh_file_name, periodicity)
    

    mpasOcean = MPAS_Ocean(mesh_directory, base_mesh_file_name, mesh_file_name, periodicity=periodicity)

    
    
    myAnalyticalGradientAtEdge = zeros(Float64, mpasOcean.nEdges, 2)
    for iEdge = 1:mpasOcean.nEdges
        myAnalyticalGradientAtEdge[iEdge,1], myAnalyticalGradientAtEdge[iEdge,2] =
                    surface_elevation_gradient(mpasOcean.lX, mpasOcean.lY, mpasOcean.xEdge[iEdge], mpasOcean.yEdge[iEdge])
    end
    myAnalyticalGradientNormalToEdge = ComputeNormalAndTangentialComponentsAtEdge(myAnalyticalGradientAtEdge, mpasOcean.angleEdge, "normal")
    
    ssh = zeros(Float64, mpasOcean.nCells)
    for iCell = 1:mpasOcean.nCells
        ssh[iCell] = surface_elevation(mpasOcean.lX, mpasOcean.lY, mpasOcean.xCell[iCell], mpasOcean.yCell[iCell])
    end
    myNumericalGradientNormalToEdge = numerical_gradient_operator(mpasOcean, ssh, periodicity)
    
    if periodicity == "NonPeriodic_x" || periodicity == "NonPeriodic_y" || periodicity == "NonPeriodic_xy"
        for iEdge = 1:mpasOcean.nEdges
            if mpasOcean.boundaryEdge[iEdge] == 1.0
                myNumericalGradientNormalToEdge[iEdge] = myAnalyticalGradientNormalToEdge[iEdge]
            end
        end
    end
    
    myGradientError = myNumericalGradientNormalToEdge - myAnalyticalGradientNormalToEdge
    
    MaxErrorNorm = norm(myGradientError, Inf)
    L2ErrorNorm = norm(myGradientError) / sqrt(mpasOcean.nEdges - mpasOcean.nNonPeriodicBoundaryEdges)
    
    @sprintf("Maximum Error Norm of the gradient is %.3g.", MaxErrorNorm)
    @sprintf("L2 Error norm of the gradient is %.3g.", L2ErrorNorm)
    
    if plotFigures
        ylabel = "Meridional Distance (km)"
        xlabel = "Zonal Distance (km)"
        
        fig, ax, cbar, _ = edgeHeatMapMesh(mpasOcean, myAnalyticalGradientNormalToEdge)
        ax.set_title("Analytical Gradient on Edges")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        vectorStreamPlotMesh(mpasOcean, myAnalyticalGradientNormalToEdge, fig=fig, ax=ax, cbar=cbar)
        display(fig)
        
#         ssh = [surface_elevation(mpasOcean.lX, mpasOcean.lY, mpasOcean.xCell[iCell], mpasOcean.yCell[iCell]) for iCell = 1:mpasOcean.nCells]
#         fig, ax, _, _ = heatMapMesh(mpasOcean, ssh)
#         ax.set_title("ssh")
#         ax.set_xlabel(xlabel)
#         ax.set_ylabel(ylabel)
#         display(fig)
        
        fig, ax, cbar, _ = edgeHeatMapMesh(mpasOcean, myNumericalGradientNormalToEdge)
        ax.set_title("Numerical Gradient on Edges")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        vectorStreamPlotMesh(mpasOcean, myNumericalGradientNormalToEdge, fig=fig, ax=ax, cbar=cbar)
        display(fig)
        
        fig, ax, cbar, _ = edgeHeatMapMesh(mpasOcean, myGradientError)
        ax.set_title("Gradient Error on Edges")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        vectorStreamPlotMesh(mpasOcean, myGradientError, fig=fig, ax=ax, cbar=cbar)
        display(fig)
    end
    
    return mpasOcean.nCells, MaxErrorNorm, L2ErrorNorm
end


# In[27]:



# <a id="gradient"></a>
# [Back to Top](#toc)

# # Divergence operator test

# In[28]:


function numerical_divergence_operator(mpasOcean, myVectorLengthNormalToEdge, periodicity)
    myDivergenceOfCell = zeros(Float64, mpasOcean.nCells)
    for iCell = 1:mpasOcean.nCells
        if !(periodicity !== "Periodic" && mpasOcean.boundaryCell[iCell] == 1.0)
            for i = 1:mpasOcean.nEdgesOnCell[iCell]
                iEdge = mpasOcean.edgesOnCell[i, iCell]
                myDivergenceOfCell[iCell] -= mpasOcean.dvEdge[iEdge] * myVectorLengthNormalToEdge[iEdge] * mpasOcean.edgeSignOnCell[iCell, i]
            end
            myDivergenceOfCell[iCell] /= mpasOcean.areaCell[iCell]
        end
    end
    return myDivergenceOfCell
end


# In[29]:


function test_numerical_divergence_operator(plotFigures,
                                    mesh_directory, base_mesh_file_name, mesh_file_name, periodicity)

    mpasOcean = MPAS_Ocean(mesh_directory, base_mesh_file_name, mesh_file_name, periodicity=periodicity)
    
    myAnalyticalSSHGradient = zeros(Float64, mpasOcean.nEdges, 2)
    for iEdge = 1:mpasOcean.nEdges
        myAnalyticalSSHGradient[iEdge,1], myAnalyticalSSHGradient[iEdge,2] =
                            surface_elevation_gradient(mpasOcean.lX, mpasOcean.lY, mpasOcean.xEdge[iEdge], mpasOcean.yEdge[iEdge])
    end
    myAnalyticalSSHGradientNormalToEdge = ComputeNormalAndTangentialComponentsAtEdge(myAnalyticalSSHGradient, mpasOcean.angleEdge, "normal")
    
    
    myNumericalSSHLaplacian = numerical_divergence_operator(mpasOcean, myAnalyticalSSHGradientNormalToEdge, periodicity)
    
    myAnalyticalSSHLaplacian = zeros(Float64, mpasOcean.nCells)
    for iCell = 1:mpasOcean.nCells
        myAnalyticalSSHLaplacian[iCell] = surface_elevation_laplacian(mpasOcean.lX, mpasOcean.lY, mpasOcean.xCell[iCell], mpasOcean.yCell[iCell])
    end
    
    
    if periodicity == "NonPeriodic_x" || periodicity == "NonPeriodic_y" || periodicity == "NonPeriodic_xy"
        for iCell = 1:mpasOcean.nCells
            if mpasOcean.boundaryCell[iCell] == 1.0
                myNumericalSSHLaplacian[iCell] = myAnalyticalSSHLaplacian[iCell]
            end
        end
    end
    
    
    mySSHLaplacianError = myNumericalSSHLaplacian - myAnalyticalSSHLaplacian
    
    MaxErrorNorm = norm(mySSHLaplacianError, Inf)
    L2ErrorNorm = norm(mySSHLaplacianError) / sqrt(mpasOcean.nCells)
    
    if plotFigures
        xlabel = "Zonal Distance (km)"
        ylabel = "Meridional Distance (km)"
        
        fig, ax, cbar, _ = heatMapMesh(mpasOcean, myAnalyticalSSHLaplacian)
        ax.set_title("Analytical Laplacian of Surface Elevation")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        vectorStreamPlotMesh(mpasOcean, myAnalyticalSSHGradientNormalToEdge, fig=fig, ax=ax, cbar=cbar)
        display(fig)
        
        
        fig, ax, cbar, _ = heatMapMesh(mpasOcean, myNumericalSSHLaplacian)
        ax.set_title("Numerical Laplacian of Surface Elevation")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        vectorStreamPlotMesh(mpasOcean, myAnalyticalSSHGradientNormalToEdge, fig=fig, ax=ax, cbar=cbar)
        display(fig)
        
        fig, ax, _, _ = heatMapMesh(mpasOcean, mySSHLaplacianError)
        ax.set_title("Numerical Error on Laplacian of Surface Elevation")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        display(fig)
    end
    
    return mpasOcean.nCells, MaxErrorNorm, L2ErrorNorm
end




# <a id="divergence"></a>
# [Back to Top](#toc)

# # Convergence Tests

# ## Tangential Velocity

# In[32]:


convergence_test("Periodic",
            CODE_ROOT * "InertiaGravityWaveMesh/ConvergenceStudyMeshes",
            "Tangential Velocity",
            test_tangential_velocity, showplots=true, savedata=true)


# ## Curl (vertex)

# In[33]:


convergence_test("Periodic",
                    CODE_ROOT * "InertiaGravityWaveMesh/ConvergenceStudyMeshes/",
                    "Numerical Curl (at Vertex)",
                    test_curl_at_vertex, savedata=true)


# ## Curl (cell center)

# In[34]:


convergence_test("Periodic",
                CODE_ROOT * "InertiaGravityWaveMesh/ConvergenceStudyMeshes/",
                "Numerical Curl (at Cell Center)",
                test_curl_at_cell_center, savedata=true)


# ## Gradient

# In[35]:


convergence_test("Periodic",
                    CODE_ROOT * "InertiaGravityWaveMesh/ConvergenceStudyMeshes/",
                    "Numerical Gradient",
                    test_numerical_gradient_operator, savedata=true)


# ## Divergence

# In[36]:


convergence_test("Periodic",
            CODE_ROOT * "InertiaGravityWaveMesh/ConvergenceStudyMeshes/",
            "Numerical Divergence",
            test_numerical_divergence_operator, savedata=true)


# In[ ]:





# In[ ]:




