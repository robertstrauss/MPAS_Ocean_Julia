include("calculate_tendencies_cuda.jl")

######## CUDA GPU versions of time-advancing  methods


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
