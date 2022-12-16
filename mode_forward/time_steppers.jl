include("calculate_tendencies.jl")

function forward_backward_step!(mpasOcean::MPAS_Ocean)
    calculate_normal_velocity_tendency!(mpasOcean)

    update_normal_velocity_by_tendency!(mpasOcean)

    calculate_ssh_tendency!(mpasOcean)

    update_ssh_by_tendency!(mpasOcean)
end

function forward_euler_step!(mpasOcean::MPAS_Ocean)
    calculate_normal_velocity_tendency!(mpasOcean)

    calculate_ssh_tendency!(mpasOcean)

    update_normal_velocity_by_tendency!(mpasOcean)

    update_ssh_by_tendency!(mpasOcean)
end

