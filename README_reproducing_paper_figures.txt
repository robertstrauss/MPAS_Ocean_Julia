Figure 1: Operator and IG wave numerical solution convergence
	plotted in output/operator_convergence/operator_convergence_plotting.ipynb
	using data from /output/operator_convergence/[...]/NonPeriodic_xy/[...].txt
		produced in /Operator_testing.ipynb
	also using plots from /InertiaGravityWaveConvergenceTest.ipynb


Table 1: 

Table 2:

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
	
