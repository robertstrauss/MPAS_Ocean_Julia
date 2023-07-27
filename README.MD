# MPAS Ocean in Julia
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7493064.svg)](https://doi.org/10.5281/zenodo.7493064)

An ocean model capable of running on irregular, non-rectilinear, TRiSK-based meshes. Inspired by MPAS-Ocean in Fortran.

**Why remake a great, working ocean model in Julia?**

Some languages are easy to develop at the cost of executing slowly, while others are lightning fast at run time but much more difficult to write. Julia is a programming language that aims to be the best of both worlds, a development and production language at the same time. To test Julia’s utility in scientific high-performance computing (HPC), we built a MPAS Shallow Water core in Julia and compared it with existing codes. We began with a simple single-processor implementation of the equations, and then ran unit tests and convergence studies for verification purpose. Subsequently, we parallelized the code by rewriting the time integrators as graphics card kernels. The GPU accelerated model achieved an amazing 500 times performance boost from the single-processor version. Finally, to make our model comparable to Fortran MPAS-Ocean, we wrote methods to divide the computational work over multiple cores of a cluster with MPI. We then performed equivalent simulations with the Julia and Fortran codes to compare the speeds, and learn how useful Julia might be for climate modeling and other HPC applications.

The main experiment done in this project are:
* Verification/validation of the mesh and equation set implementation.
* Strong and weak scaling tests run on high-performance clusters, using MPI to distribute simulation across many nodes
* Performance benchmarking of the simulation using the graphics card (GPU), compared to distributed simulation with MPI.

Currently only includes gravity, coriolis terms (no non-linear terms).


# How to use and reproduce results:

## Requirements
1. Install the Julia language (we tested with version 1.8.3): https://julialang.org/downloads/
2. Many of our experiements are done within Jupyter Notebooks, for convenience of viewing graphical output and code together interacively. Install jupyter notebook/lab: https://jupyter.org/install
3. To run the distributed (MPI) or graphics card (GPU) simulations will require an MPI-compatible cluster or NVIDIA GPU respectively, and corresponding libraries and drivers. We tested with MPICH_jll 4.0.2, CUDA toolkit 11.7 and NVIDIA driver 525.105.17. Data from running the simulations on NERSC's Perlmutter cluster (https://docs.nersc.gov/systems/perlmutter/architecture/) is included. Serial versions of the simulation can be run with just a CPU.

## Preperation
4. Clone the repository and open it <br>
    $ `git clone git@github.com:/robertstrauss/MPAS_Ocean_Julia` or `git clone https://github.com/robertstrauss/MPAS_Ocean_Julia.git` <br>
    $ `cd MPAS_Ocean_Julia` <br>
5. Install the required Julia packages (may take some time to download) <br>
    $ `julia` <br>
    julia> `] activate .` <br>
    pkg> `instantiate` <br>
6. Download the mesh files from https://zenodo.org/record/7411962, (extract the zip file), and place the `MPAS_Ocean_Shallow_Water_Meshes_Julia_Paper/InertiaGravityWaveMesh/ConvergenceStudyMeshes/` directory at `./MPAS_O_Shallow_Water/ConvergenceStudyMeshes/InertiaGravityWave`, and the `MPAS_Ocean_Shallow_Water_Meshes_Julia_Paper/CoastalKelvinWaveMesh/ConvergenceStudyMeshes` directory at `./MPAS_O_Shallow_Water/ConvergenceStudyMeshes/CoastalKelvinWave` relative to this repository. Also place the `InertiaGravityWaveMesh/` director in the root of this repository, and `CoastalKelvinWaveMesh/` in the root of this repository.<br>
    $ `wget -O MPAS_Ocean_Shallow_Water_Meshes.zip https://zenodo.org/record/7411962/files/MPAS_Ocean_Shallow_Water_Meshes.zip?download=1` <br>
    $ `unzip MPAS_Ocean_Shallow_Water_Meshes.zip` <br>
    $ `mkdir -p ./MPAS_O_Shallow_Water/ConvergenceStudyMeshes/` <br>
    $ `cp -R MPAS_Ocean_Shallow_Water_Meshes/MPAS_Ocean_Shallow_Water_Meshes_Julia_Paper/InertiaGravityWaveMesh/ConvergenceStudyMeshes/ ./MPAS_O_Shallow_Water/ConvergenceStudyMeshes/InertiaGravityWave/` <br>
    $ `cp -r MPAS_Ocean_Shallow_Water_Meshes/InertiaGravityWaveMesh/ .` <br>
    $ `cp -R MPAS_Ocean_Shallow_Water_Meshes/MPAS_Ocean_Shallow_Water_Meshes_Julia_Paper/CoastalKelvinWaveMesh/ConvergenceStudyMeshes/ ./MPAS_O_Shallow_Water/ConvergenceStudyMeshes/CoastalKelvinWave` <br>
    $ `cp -r MPAS_Ocean_Shallow_Water_Meshes/CoastalKelvinWaveMesh/ .` <br>
7. For rerunning the distributed simulation tests with MPI, install the julia version of the tool for starting a distributed script across nodes: <br>
    julia> `using MPI` <br>
    julia> `MPI.install_mpiexecjl()` <br>
    julia> `exit()` <br>
    a. (Optional) add the tool to your path for easier use: <br>
    $ `ln -s ~/.julia/bin/mpiexecjl <some place on your $PATH>` <br>
    
## Reproducing figures from the paper

To reproduce the figures in the paper:
 * Figure 1 (Operator and Inertia-Gravity-Wave Convergence): run the Julia files `Operator_testing.jl` to run the tests and generate the data, and `followed by `operator_convergence_plotting.jl` to create the plots from this data at `./output/operator_convergence/<operator>/Periodic/<figure>.pdf`. Run the Julia script `InertiaGravityWaveConvergenceTest.jl` to generate the numerical solution time-horizon seas-surface height at `./output/simulation_convergence/inertiagravitywave/Periodic/CPU/` and convergence plot. For more information, visuals, and intermediate tests, use the Jupyter notebook versions of these scripts: `./Operator_testing.ipynb`, `./operator_convergence_plotting.ipynb`, and `./InertiaGravityWaveConvergenceTest.ipynb`. (A Jupyter Notebook with the proper Julia kernel can be started by running `using IJulia; IJulia.Notebook()` in the Julia shell `julia --project=.`) <br>
    $ `julia --project=. ./Operator_testing.jl` <br>
    $ `julia --project=. ./operator_convergence_plotting.jl` <br>
    $ `julia --project=. ./InertiaGravityWaveConvergenceTest.jl` <br>
 * Figures 2 (Scaling on One Node, MPI and GPU), 3 (Strong Scaling), 4 (Weak Scaling), and 5 (Time Proportion): on a cluster with at least 128 nodes with 64 processes per node, use the script `./run_scaling_16x_to_512x.sh` after modifying the slurm options specified appropriately for your cluster to execute the performance scaling tests on each mesh resolution from 16x16 to 512x512. This will save its results in `./output/kelvinwave/resolution<mesh size>/procs<maximum number of processors>/steps10/nvlevels100/`. Also execute the julia script `GPU_performance.jl` on a node with an NVIDIA graphics card to do performance tests on the GPU. Then, execute the notebook `./scalingplots.ipynb` or the julia script `scalingplots.jl` to create the figures in the paper at `./plots/<type>/<figure>.pdf`. Alternatively, use the jupyter notebooks `GPU_Performance.ipynb` and `scalingplots.ipynb` for more information, intermediate tests, and visuals. <br>
    $ `vim ./run_scaling_16x_to_512x.sh` <br>
    $ `bash ./run_scaling_16x_to_512x.sh` <br>
    $ `salloc -C gpu_count:1 -t 20:00`
    (gpu node)$ `julia --project=. ./GPU_Performance.jl` <br>
    (gpu node)$ `exit`
    $ `julia --project=. ./scalingplots.jl` <br>
 * Tables 1 & 2: execute `./serial_julia_performance.jl` with julia to produce the optimized serial timing data. Then, download the unoptimized version of the codebase from https://github.com/robertstrauss/MPAS_Ocean_Julia/tree/unoptimized or MPAS_Ocean_Julia-unopt.zip from https://doi.org/10.5281/zenodo.7493065 , and execute the julia script `./serial_julia_performance.jl` in the directory of the unoptimized codebase as well. The results will be saved in text files at `./output/serialCPU_timing/coastal_kelvinwave/unoptimized/steps_10/resolution_<mesh size>/` in the unoptimized directory and `./output/serialCPU_timing/coastal_kelvinwave/steps_10/resolution_<mesh size>/` in the main/optimized directory. <br>
    $ `julia --project=. ./serial_julia_performance.jl <mesh size x> <number of samples>` <br>
    (if cloned with git)$ `git checkout unoptimized` <br>
    $ `julia --project=. ./serial_julia_performance.jl <mesh size x> <number of samples>` <br>
    


# How to run a simulation with MPI (distributed) scaling test
Data from tests run on NERSC's cori-haswell and perlmutter are included in the `./output/` directory.
To re-run the simulations on your cluster and create new data:
* Use `./scaling_test.jl` and specify the tests to do: <br>
    $ `~/.julia/bin/mpiexecjl --project=. -n <nprocs> julia ./scaling_test/scaling_test.jl <cells> <samples> <proccounts>` <br>
    Where \<nprocs\> = the number of processors needed for the trial distributed over the most processors (should be power of 2), <br>
        \<cells\> = the width of the mesh to use for all trials (128, 256, 512 used in paper), <br>
        \<samples\> = the number of times to re-run a simulation for each trial (reccomended at least 2, the first run is often an outlier), <br>
        \<proccounts\> = a list of how many processor should be used for each trial (use powers of 2, seperated by commas, e.g. 2,8,16,32,256 would run 5 trials, the first distributed over 2 procs, then 8, then 16, etc....) Leave blank to run all powers of 2 up to the number of allocated processors (\<nprocs\>). <br>
    - you may need to allocate this many processors on the cluster you use to perform this.
* The results will be saved in `./output/kelvinwave/resolution<cells>x<cells>/procs<nprocs>/steps10/nvlevels100/` <br>

# Verifying/validating the model
To verify the implementation of the mesh, see `./Operator_testing.ipynb`.
To validate the implementation of the shallow water equations, see `InertiaGravityWaveConvergenceTest.ipynb`.
