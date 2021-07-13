# architecture scheme inspired by Oceananigans.jl
# https://github.com/CliMA/Oceananigans.jl
# having a CPU or GPU struct allows us to overload functions
# based on what architecture or device the model is on
# this makes it very easy to quickly change the device a model uses
# without changing all the time stepping methods or anything else

# abstract type AbstractArchitecture end

# abstract type AbstractCPUArchitecture <: AbstractArchitecture end

# abstract type AbstractGPUArchitecture <: AbstractArchitecture end

# struct CPU <: AbstractCPUArchitecture end

# struct GPU <: AbstractGPUArchitecture end

# device(::AbstractCPUArchitecture) = KernelAbstractions.CPU()
# device(::AbstractGPUArchitecture) = CUDAKernels.CUDADevice()

# array_type(::CPU) = Array
# array_type(::GPU) = CUDA.CuArray


# abstract type RelocatableSimulation{A<:AbstractArchitecture}
#     # Relocatable Simulations (RS) are models that can be moved between the GPU or CPU
# end


# function relocated(rs::RelocatableSimulation, architecture)
    
#     fieldvalues = []
    
#     for field in fieldnames(typeof(rs))
#         if typeof(getfield(rs, field)) <: AbstractArray
#             val = array_type(architecture)(getfield(rs, field))
#         else
#             val = getfield(rs, field)
#         end
        
#         append!(fieldvalues, [val])
#     end
    
#     return eval(typeof(rs).name.name){typeof(architecture)}(fieldvalues...)
# end