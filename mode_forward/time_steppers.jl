include("calculate_tendencies.jl")
include("../mode_init/MPAS_OceanHalos.jl")

function forward_backward_step!(mpasOceanHalos::MPAS_OceanHalos)
    for mpasOcean in mpasOceanHalos.haloChunks
	forward_backward_step!(mpasOcean)
    end
end

function forward_backward_step!(mpasOcean::MPAS_Ocean)
    calculate_normal_velocity_tendency!(mpasOcean)
    
    update_normal_velocity_by_tendency!(mpasOcean)
    
    calculate_ssh_tendency!(mpasOcean)
    
    update_ssh_by_tendency!(mpasOcean)
end

function forward_backward_step_cuda!(mpasOcean::MPAS_Ocean)
    calculate_normal_velocity_tendency_cuda!(mpasOcean)
    
    update_normal_velocity_by_tendency_cuda!(mpasOcean)
    
    calculate_ssh_tendency_cuda!(mpasOcean)
    
    update_ssh_by_tendency_cuda!(mpasOcean)
end







function forward_euler_step!(mpasOcean::MPAS_Ocean)
    calculate_normal_velocity_tendency!(mpasOcean)
    
    calculate_ssh_tendency!(mpasOcean)
    
    update_normal_velocity_by_tendency!(mpasOcean)

    update_ssh_by_tendency!(mpasOcean)
end

function forward_euler_step_cuda!(mpasOcean::MPAS_Ocean)
    calculate_normal_velocity_tendency_cuda!(mpasOcean)
    
    calculate_ssh_tendency_cuda!(mpasOcean)
    
    update_normal_velocity_by_tendency_cuda!(mpasOcean)

    update_ssh_by_tendency_cuda!(mpasOcean)
end
