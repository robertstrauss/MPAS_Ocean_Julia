
mutable struct Namelist

    config_bottom_depth_parameter

    config_mesh_type
    config_mesh_problem

    config_problem_type_Geophysical_Wave
    config_problem_type_Equatorial_Wave

    config_problem_is_linear
    config_periodicity
    config_time_integrator

    config_Forward_Backward_with_RK2_Feedback_parameter_beta
    config_Forward_Backward_with_RK2_Feedback_parameter_epsilon
    config_LF_TR_and_LF_AM3_with_FB_Feedback_Type

    config_LF_TR_and_LF_AM3_with_FB_Feedback_parameter_beta
    config_LF_TR_and_LF_AM3_with_FB_Feedback_parameter_gamma
    config_LF_TR_and_LF_AM3_with_FB_Feedback_parameter_epsilon

    config_Generalized_FB_with_AB2_AM3_Step_Type

    config_Generalized_FB_with_AB2_AM3_Step_parameter_beta
    config_Generalized_FB_with_AB2_AM3_Step_parameter_gamma
    config_Generalized_FB_with_AB2_AM3_Step_parameter_epsilon

    config_Generalized_FB_with_AB3_AM4_Step_Type

    config_Generalized_FB_with_AB3_AM4_Step_parameter_beta
    config_Generalized_FB_with_AB3_AM4_Step_parameter_gamma
    config_Generalized_FB_with_AB3_AM4_Step_parameter_epsilon
    config_Generalized_FB_with_AB3_AM4_Step_parameter_delta



    config_bottom_depth_parametermyNamelist
    config_bottom_slope
    config_Coriolis_parameter
    config_meridional_gradient_of_Coriolis_parameter
    config_gravity
    config_mean_depth
    config_surface_elevation_amplitude
    config_problem_type

    config_thickness_flux_type
    config_use_wetting_drying
    # Derived Parameters
    config_phase_speed_of_coastal_Kelvin_wave
    config_radius_of_deformation
    config_equatorial_radius_of_deformation
    config_dt

    config_linearity_prefactor

    config_zonal_diffusivity

    config_meridional_diffusivity
    config_viscous_Burgers_zonal_velocity_left
    config_viscous_Burgers_zonal_velocity_right
    config_viscous_Burgers_shock_speed
    config_viscous_Burgers_zonal_offset




    



    function Namelist(;mesh_type="uniform",problem_type="default",problem_is_linear=true,
                 periodicity="Periodic",time_integrator="Forward_Backward",
                 LF_TR_and_LF_AM3_with_FB_Feedback_Type="FourthOrderAccurate_MaximumStabilityRange",
                 Generalized_FB_with_AB2_AM3_Step_Type="FourthOrderAccurate",
                 Generalized_FB_with_AB3_AM4_Step_Type="FourthOrderAccurate_MaximumStabilityRange")
        myNamelist = new()
        myNamelist.config_mesh_type = mesh_type
        myNamelist.config_problem_type = problem_type
        if (problem_type == "Coastal_Kelvin_Wave" || problem_type == "Inertia_Gravity_Wave"
            || problem_type == "Inertia_Gravity_Waves" || problem_type == "Planetary_Rossby_Wave"
            || problem_type == "Topographic_Rossby_Wave" || problem_type == "Equatorial_Kelvin_Wave"
            || problem_type == "Equatorial_Yanai_Wave" || problem_type == "Equatorial_Rossby_Wave"
            || problem_type == "Equatorial_Inertia_Gravity_Wave")
            myNamelist.config_problem_type_Geophysical_Wave = true
        else
            myNamelist.config_problem_type_Geophysical_Wave = false
        end
        if (problem_type == "Equatorial_Kelvin_Wave" || problem_type == "Equatorial_Yanai_Wave"
            || problem_type == "Equatorial_Rossby_Wave" || problem_type == "Equatorial_Inertia_Gravity_Wave")
            myNamelist.config_problem_type_Equatorial_Wave = true
        else
            myNamelist.config_problem_type_Equatorial_Wave = false
        end
        myNamelist.config_problem_is_linear = problem_is_linear
        myNamelist.config_periodicity = periodicity
        myNamelist.config_time_integrator = time_integrator
        if time_integrator == "Forward_Backward_with_RK2_Feedback"
            myNamelist.config_Forward_Backward_with_RK2_Feedback_parameter_beta = 1.0/3.0
            myNamelist.config_Forward_Backward_with_RK2_Feedback_parameter_epsilon = 2.0/3.0
        elseif time_integrator == "LF_TR_and_LF_AM3_with_FB_Feedback"
            myNamelist.config_LF_TR_and_LF_AM3_with_FB_Feedback_Type = LF_TR_and_LF_AM3_with_FB_Feedback_Type
            if LF_TR_and_LF_AM3_with_FB_Feedback_Type == "SecondOrderAccurate_LF_TR"
                beta = 0.0
                gamma = 0.0
                epsilon = 0.0
            elseif LF_TR_and_LF_AM3_with_FB_Feedback_Type == "ThirdOrderAccurate_LF_AM3"
                beta = 0.0
                gamma = 1.0/12.0
                epsilon = 0.0
            elseif LF_TR_and_LF_AM3_with_FB_Feedback_Type == "ThirdOrderAccurate_MaximumStabilityRange"
                beta = 0.126
                gamma = 1.0/12.0
                epsilon = 0.83
            elseif LF_TR_and_LF_AM3_with_FB_Feedback_Type == "FourthOrderAccurate_MinimumTruncationError"
                beta = 17.0/120.0
                gamma = 1.0/12.0
                epsilon = 11.0/20.0
            elseif LF_TR_and_LF_AM3_with_FB_Feedback_Type == "FourthOrderAccurate_MaximumStabilityRange"
                epsilon = 0.7166
                beta = 7.0/30.0 - epsilon/6.0
                gamma = 1.0/12.0
            end
            myNamelist.config_LF_TR_and_LF_AM3_with_FB_Feedback_parameter_beta = beta
            myNamelist.config_LF_TR_and_LF_AM3_with_FB_Feedback_parameter_gamma = gamma
            myNamelist.config_LF_TR_and_LF_AM3_with_FB_Feedback_parameter_epsilon = epsilon
        elseif time_integrator == "Generalized_FB_with_AB2_AM3_Step"
            myNamelist.config_Generalized_FB_with_AB2_AM3_Step_Type = Generalized_FB_with_AB2_AM3_Step_Type
            if Generalized_FB_with_AB2_AM3_Step_Type == "ThirdOrderAccurate_WideStabilityRange"
                beta = 0.0
            elseif (Generalized_FB_with_AB2_AM3_Step_Type
                  == "ThirdOrderAccurate_WeakAsymptoticInstabilityOfPhysicalModes")
                beta = 0.5
            elseif Generalized_FB_with_AB2_AM3_Step_Type == "FourthOrderAccurate"
                useSymPyToDetermineBeta = false
                # Note that if useSymPyToDetermineBeta is specified as true, every time beta is used, the SymPy
                # polynomial equation solver will be executed, resulting in immense slowdown of the code.
                if useSymPyToDetermineBeta
                    symbolic_beta = sp.Symbol("beta")
                    beta_roots = sp.solve(-symbolic_beta^3.0 - symbolic_beta/12.0 + 1.0/12.0, symbolic_beta)
                    beta = beta_roots[1]
                else
                    beta = 0.373707625197906
                end
            end
            gamma = beta - 2.0*beta^2.0 - 1.0/6.0
            epsilon = beta^2.0 + 1.0/12.0
            myNamelist.config_Generalized_FB_with_AB2_AM3_Step_parameter_beta = beta
            myNamelist.config_Generalized_FB_with_AB2_AM3_Step_parameter_gamma = gamma
            myNamelist.config_Generalized_FB_with_AB2_AM3_Step_parameter_epsilon = epsilon
        end
        if time_integrator == "Generalized_FB_with_AB3_AM4_Step"
            myNamelist.config_Generalized_FB_with_AB3_AM4_Step_Type = Generalized_FB_with_AB3_AM4_Step_Type
            if Generalized_FB_with_AB3_AM4_Step_Type == "SecondOrderAccurate_OptimumChoice_ROMS"
                beta = 0.281105
                gamma = 0.088
                epsilon = 0.013
            elseif Generalized_FB_with_AB3_AM4_Step_Type == "ThirdOrderAccurate_AB3_AM4"
                beta = 5.0/12.0
                gamma = -1.0/12.0
                epsilon = 0.0
            elseif Generalized_FB_with_AB3_AM4_Step_Type == "ThirdOrderAccurate_MaximumStabilityRange"
                beta = 0.232
                epsilon = 0.00525
                gamma = 1.0/3.0 - beta - 3.0*epsilon
            elseif Generalized_FB_with_AB3_AM4_Step_Type == "ThirdOrderAccurate_OptimumChoice"
                beta = 0.21
                epsilon = 0.0115
                gamma = 1.0/3.0 - beta - 3.0*epsilon
            elseif Generalized_FB_with_AB3_AM4_Step_Type == "FourthOrderAccurate_MaximumStabilityRange"
                epsilon = 0.083
                gamma = 0.25 - 2.0*epsilon
                beta = 1.0/12.0 - epsilon
            end
            delta = 0.5 + gamma + 2.0*epsilon
            myNamelist.config_Generalized_FB_with_AB3_AM4_Step_parameter_beta = beta
            myNamelist.config_Generalized_FB_with_AB3_AM4_Step_parameter_gamma = gamma
            myNamelist.config_Generalized_FB_with_AB3_AM4_Step_parameter_epsilon = epsilon
            myNamelist.config_Generalized_FB_with_AB3_AM4_Step_parameter_delta = delta
        end
        myNamelist.config_bottom_depth_parameter = 1000.0
        myNamelist.config_bottom_slope = 0.0
        myNamelist.config_Coriolis_parameter = 10.0^(-4.0)
        myNamelist.config_meridional_gradient_of_Coriolis_parameter = 2.0*10.0^(-11.0)
        myNamelist.config_gravity = 10.0
        myNamelist.config_mean_depth = 1000.0
        if (problem_type == "default" || problem_type == "Coastal_Kelvin_Wave"
            || myNamelist.config_problem_type_Equatorial_Wave)
            myNamelist.config_surface_elevation_amplitude = 0.001
        elseif (problem_type == "Inertia_Gravity_Wave" || problem_type == "Planetary_Rossby_Wave"
              || problem_type == "Topographic_Rossby_Wave" || problem_type == "Diffusion_Equation")
            myNamelist.config_surface_elevation_amplitude = 1.0
        elseif myNamelist.config_problem_type == "Barotropic_Tide"
            myNamelist.config_surface_elevation_amplitude = 2.0
        elseif myNamelist.config_problem_type == "Viscous_Burgers_Equation"
            myNamelist.config_surface_elevation_amplitude = 0.0
        end
        myNamelist.config_thickness_flux_type = "centered"
        myNamelist.config_use_wetting_drying = false
        # Derived Parameters
        myNamelist.config_phase_speed_of_coastal_Kelvin_wave = (
        sqrt(myNamelist.config_gravity*myNamelist.config_mean_depth))
        myNamelist.config_radius_of_deformation = (
        myNamelist.config_phase_speed_of_coastal_Kelvin_wave/myNamelist.config_Coriolis_parameter)
        myNamelist.config_equatorial_radius_of_deformation = (
        sqrt(myNamelist.config_phase_speed_of_coastal_Kelvin_wave
                /myNamelist.config_meridional_gradient_of_Coriolis_parameter))
        if problem_type == "default" || problem_type == "Coastal_Kelvin_Wave"
            myNamelist.config_dt = 180.0
        elseif problem_type == "Inertia_Gravity_Wave"
            myNamelist.config_dt = 96.0
        elseif problem_type == "Planetary_Rossby_Wave" || problem_type == "Topographic_Rossby_Wave"
            myNamelist.config_dt = 195000.0
        elseif problem_type == "Equatorial_Kelvin_Wave"
            myNamelist.config_dt = 750.00
        elseif problem_type == "Equatorial_Yanai_Wave"
            myNamelist.config_dt = 390.00
        elseif problem_type == "Equatorial_Rossby_Wave"
            myNamelist.config_dt = 2700.00
        elseif problem_type == "Equatorial_Inertia_Gravity_Wave"
            myNamelist.config_dt = 420.00
        elseif problem_type == "Barotropic_Tide"
            myNamelist.config_dt = 10.00
        elseif problem_type == "Diffusion_Equation"
            myNamelist.config_dt = 2260.0
        elseif problem_type == "Viscous_Burgers_Equation"
            myNamelist.config_dt = 2100.00
        end
        if problem_is_linear
            myNamelist.config_linearity_prefactor = 0.0
        else
            myNamelist.config_linearity_prefactor = 1.0
        end
        if problem_type == "Viscous_Burgers_Equation"
            myNamelist.config_zonal_diffusivity = 5000.0
        else
            myNamelist.config_zonal_diffusivity = 500.0
        end
        myNamelist.config_meridional_diffusivity = 500.0
        myNamelist.config_viscous_Burgers_zonal_velocity_left = 1.0
        myNamelist.config_viscous_Burgers_zonal_velocity_right = 0.0
        myNamelist.config_viscous_Burgers_shock_speed = (
        0.5*(myNamelist.config_viscous_Burgers_zonal_velocity_left
             + myNamelist.config_viscous_Burgers_zonal_velocity_right))
        myNamelist.config_viscous_Burgers_zonal_offset = 0.0

        return myNamelist
    end
end
