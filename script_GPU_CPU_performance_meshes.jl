# coding: utf-8

# In[1]:


CODE_ROOT = pwd() * "/"

include(CODE_ROOT * "mode_forward/time_steppers.jl")
include(CODE_ROOT * "mode_init/MPAS_Ocean.jl")
include(CODE_ROOT * "visualization.jl")
include(CODE_ROOT * "mode_init/exactsolutions.jl")
include(CODE_ROOT * "mode_init/MPAS_Ocean_CUDA.jl")


# In[2]:


using BenchmarkTools
using DelimitedFiles
using Dates


# In[23]:


function kelvinWaveCUDAGenerator(mpasOcean, lateralProfile)
    meanCoriolisParameterf = sum(mpasOcean.fEdge) / length(mpasOcean.fEdge)
    meanFluidThicknessH = sum(mpasOcean.bottomDepth)/length(mpasOcean.bottomDepth)
    c = sqrt(mpasOcean.gravity*meanFluidThicknessH)
    rossbyRadiusR = c/meanCoriolisParameterf
    
    # writing the exact solution to the boundary edges
    function boundaryConditionCUDA_kernel!(nEdges, normalVelocityCurrent, boundaryEdge, angleEdge, yEdge, xEdge, t)
        iEdge = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
        if iEdge <= nEdges
            if boundaryEdge[iEdge] == 1
                normalVelocityCurrent[iEdge,:] .= CUDA.sin(angleEdge[iEdge]) * c * lateralProfile(yEdge[iEdge] + c*t) * CUDA.exp(-xEdge[iEdge]/rossbyRadiusR)
            end
        end
        return nothing
    end
    function boundaryConditionCUDA!(mpasOcean::MPAS_Ocean_CUDA, t::Float64)
        CUDA.@cuda blocks=cld(mpasOcean.nEdges, 1024) threads=1024 maxregs=64 boundaryConditionCUDA_kernel!(mpasOcean.nEdges,
                                                                                        mpasOcean.normalVelocityCurrent,
                                                                                        mpasOcean.boundaryEdge,
                                                                                        mpasOcean.angleEdge,
                                                                                        mpasOcean.yEdge,
                                                                                        mpasOcean.xEdge,
                                                                                        t)
        return nothing
    end

    return boundaryConditionCUDA!
end


# In[83]:


function gpucpuperftest(nCellsX, nSamples=10, nSteps=10, device="CPU", showplots=false, writedata=false; nvlevels=100)
    fpath = CODE_ROOT * "output/asrock/serial$(device)_timing/coastal_kelvinwave/steps_$nSteps/resolution_$(nCellsX)x$(nCellsX)/nvlevels_$(nvlevels)/"
    fname = "$fpath$(Dates.now()).txt"
    if writedata
        mkpath(fpath)
        println("output file: $fname")
    end
    
    mpasMesh = MPAS_Ocean(CODE_ROOT * "ConvergenceStudyMeshes",
        "culled_mesh_$(nCellsX)x$(nCellsX).nc", "mesh_$(nCellsX)x$(nCellsX).nc", periodicity="NonPeriodic_x", nvlevels=nvlevels)
    mpasOcean = mpasMesh
    device_step_nv! = function(mpasOcean)
        calculate_normal_velocity_tendency!(mpasOcean)
        update_normal_velocity_by_tendency!(mpasOcean)
    end
    device_step_ssh! = function(mpasOcean)
        calculate_ssh_tendency!(mpasOcean)
        update_ssh_by_tendency!(mpasOcean)
    end

    @inline function lateralProfile(y)
        return 1e-3 * exp( - (y-3e6)^2 / 5e5^2 )
    end
    
    kelvinWaveExactNormalVelocity, kelvinWaveExactSSH, kelvinWaveExactSolution!, device_boundary_condition! = kelvinWaveGenerator(mpasOcean, lateralProfile)
    
    
    
    exacttime = 0.0
    kelvinWaveExactSolution!(mpasOcean, exacttime)
    
    if showplots
        fig, ax, _ = heatMapMesh(mpasMesh, mpasOcean.sshCurrent)
        ax.set_title("Initial Sea Surface Height")
        display(fig)
    end
    
    if "GPU" == device
        mpasOcean = MPAS_Ocean_CUDA(mpasOcean)
        device_boundary_condition! = kelvinWaveCUDAGenerator(mpasOcean, lateralProfile)
        device_step_nv! = function(mpasOcean)
            calculate_normal_velocity_tendency_cuda!(mpasOcean)
            update_normal_velocity_by_tendency_cuda!(mpasOcean)
        end
        device_step_ssh! = function(mpasOcean)
            calculate_ssh_tendency_cuda!(mpasOcean)
            update_ssh_by_tendency_cuda!(mpasOcean)
        end
    end
    
    sampletimes = zeros(Float64, nSamples)
    
    for s in 1:nSamples
        sampletimes[s] += @elapsed begin for j in 1:nSteps
                device_step_nv!(mpasOcean)
                device_boundary_condition!(mpasOcean, exacttime)
                device_step_ssh!(mpasOcean)
                exacttime += mpasOcean.dt
            end
        end
    end
    
    
    println("calculating error")
    exactSSH = zeros(Float64, mpasOcean.nCells)
    for iCell in 1:mpasOcean.nCells
        exactSSH[iCell] = kelvinWaveExactSSH(mpasMesh, iCell, exacttime)
    end
    
    error =  sum( ( Array(mpasOcean.sshCurrent) .- exactSSH ) .^ 2 )
    
    if showplots
        fig, ax, _ = heatMapMesh(mpasMesh, Array(mpasOcean.sshCurrent))
        ax.set_title("Final Numerical Sea Surface Height")
        display(fig)
        
        fig, ax, _ = heatMapMesh(mpasMesh, exactSSH)
        ax.set_title("Final Exact Sea Surface Height")
        display(fig)
    end
    
    if writedata
        open(fname, "w") do io
            writedlm(io, sampletimes)
        end
    end
    
    return sampletimes, error, mpasOcean, mpasMesh
end


# In[ ]:


resolutions = [64, 144, 324]


# In[87]:


nCellsX = resolutions[1]
wcgpu, errgpu, mpasOceanCuda, mpasMesh = gpucpuperftest(nCellsX, 10, 20, "GPU", true, true, nvlevels=100)
println("execution time: $wcgpu")
println("error: $errgpu")


# In[88]:



wccpu, errcpu, mpasOceanCuda, mpasMesh = gpucpuperftest(nCellsX, 10, 20, "CPU", true, true, nvlevels=100)
println("execution time: $wccpu")
println("error: $errcpu")


# In[90]:


avgcpu = sum(wccpu)/length(wccpu)
avggpu = sum(wcgpu)/length(wcgpu)
avgcpu, avggpu, " speedup: ", avgcpu/avggpu


# In[84]:


nCellsX = resolutions[1]
wcgpu, errgpu, mpasOceanCuda, mpasMesh = gpucpuperftest(nCellsX, 10, 20, "GPU", true, true, nvlevels=2)
println("execution time: $wcgpu")
println("error: $errgpu")


# In[ ]:


nCellsX = resolutions[1]
wccpu, errcpu, mpasOceanCuda, mpasMesh = gpucpuperftest(nCellsX, 10, 20, "CPU", true, true, nvlevels=2)
println("execution time: $wccpu")
println("error: $errcpu")


# In[4]:


for nCellsX in resolutions
    wc, err = gpucpuperftest(nCellsX, 2, 10, "CPU", true, true)
    println("wallclock $wc")
    println("error $err")
end


# In[5]:


for nCellsX in resolutions
    wc, err = gpucpuperftest(nCellsX, 2, 10, "GPU", true, true)
    println("wallclock $wc")
    println("error $err")
end


# In[11]:


nCellsX = resolutions[end]
wc, err, mpasOceanCuda = gpucpuperftest(nCellsX, 2, 10, "GPU", true, true)
mpasOcean = MPAS_Ocean(CODE_ROOT * "MPAS_O_Shallow_Water/MPAS_O_Shallow_Water_Mesh_Generation/CoastalKelvinWaveMesh/ConvergenceStudyMeshes",
        "culled_mesh_$nCellsX.nc", "mesh_$nCellsX.nc", periodicity="NonPeriodic_x", nvlevels=100)
0


# In[13]:


fig, ax, _ = heatMapMesh(mpasOcean, Array(mpasOceanCuda.sshCurrent))
display(fig)


# the data copying step is actually 2 - 5 times slower than it is on the CPU! Memory is the bottleneck of the GPU
