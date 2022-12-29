Figure 1: Operator and IG wave numerical solution convergence
    plotted in output/operator_convergence/operator_convergence_plotting.ipynb
    using data from /output/operator_convergence/[...]/NonPeriodic_xy/[...].txt
        produced in /Operator_testing.ipynb
            using mesh files MPAS_O_Shallow_Water_Mesh_Generation/ConvergenceStudyMeshes
                available at MPAS_O_Shallow_Water_Meshes/ConvergenceStudyMeshes from
                    zenodo: https://zenodo.org/record/7419817#.Y63p4C-B1pQ
    also using plots from /InertiaGravityWaveConvergenceTest.ipynb
        using mesh files from InertiaGravityWaveMesh/
            available at: zenodo: https://zenodo.org/record/7419817#.Y63p4C-B1pQgith

Table 1 and Table 2: 
    Average of data from:
        Python: Python Rotating Shallow Water Verification Suite by Siddartha Bishnu, available at:
            github: https://github.com/siddharthabishnu/Rotating_Shallow_Water_Verification_Suite.git
            zenodo: https://doi.org/10.5281/zenodo.7425628
        Julia:
            data saved at 
                output/serialCPU_timing/coastal_kelvinwave/unoptimized/steps_$nSteps/resolution_$(nCellsX)x$(nCellsX)/
                output/serialCPU_timing/coastal_kelvinwave/steps_$nSteps/resolution_$(nCellsX)x$(nCellsX)/
                output/asrock/serial$(device)_timing/coastal_kelvinwave/steps_$nSteps/resolution_$(nCellsX)x$(nCellsX)/nvlevels_$(nvlevels)/
            made with:
                serial_julia_performance.jl\t
                cpu-unoptimized:
                    github: https://github.com/robertstrauss/MPAS_Ocean_Julia/tree/unoptimized
                    zenodo: https://doi.org/10.5281/zenodo.7493065 (MPAS_Ocean_Julia-1.0.zip)
                cpu-optimized and Julia gpu:
                    GPU_CPU_performance_comparison_meshes.ipynb at
                    github: https://github.com/robertstrauss/MPAS_Ocean_Julia/tree/main
                    zenodo: https://doi.org/10.5281/zenodo.7493065 (MPAS_Ocean_Julia-2.0.zip)
            
            using mesh files from MPAS_O_Shallow_Water_Mesh_Generation/CoastalKelvinWaveMesh/ConvergenceStudyMeshes
                    available at MPAS_O_Shallow_Water_Meshes/CoastalKelvinWaveMesh/ConvergenceStudyMeshes from
                        zenodo: https://zenodo.org/record/7419817#.Y63p4C-B1pQ 

Figure 2:
    Produced by M. Petersen


Figure 3: Simulation speed and Occupancy
    plotted in /CUDAOccupancyTest.ipynb
    data saved at /output/cuda_occupancy/[...].txt

Figure 4: Scaling Test comparing julia and fortran on three mesh resolutions
and Figure 5: Proportion computation vs communication
    plotted in output/kelvinwave/performanceplots.ipynb
    using data from
        Fortran: /output/kelvinwave/fortranperformance/resolution<X>x<X>/cori-haswell_<p>_[...].txt
            Produced by M. Petersen, using MPAS on cori-haswell
            nonlinear terms disabled, other modifications for fair comparison
        Julia: /output/kelvinwave/resolution<X>x<X>/procs2048/steps10/nvlevels100/[...].txt
            Produced by running /scaling_test/scaling_test.jl on cori-haswell:
                sbatch [...] mpiexecjl -n <p> julia --project ./scaling_test/scaling_test.jl <X> <s>
                sbatch -o sbatch_logs/%A.out -e sbatch_logs/%A_%t.err -A e3sm -C haswell --distribution=plane=32 --ntasks-per-node=32 --qos=regular -t 20:00 -n 2048 mpiexecjl -n 2048 julia --project ./scaling_test/scaling_test.jl 128 12
    
