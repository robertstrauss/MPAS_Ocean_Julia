# MPAS Ocean in Julia
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7493065.svg)](https://doi.org/10.5281/zenodo.7493065)

An ocean model capable of running on irregular, non-rectilinear, TRiSK-based meshes. Inspired by MPAS-Ocean in Fortran.

**Why remake a great, working ocean model in Julia?**

Some languages are easy to develop at the cost of executing slowly, while others are lightning fast at run time but much more difficult to write. Julia is a programming language that aims to be the best of both worlds, a development and production language at the same time. To test Julia’s utility in scientific high-performance computing (HPC), we built a MPAS Shallow Water core in Julia and compared it with existing codes. We began with a simple single-processor implementation of the equations, and then ran unit tests and convergence studies for verification purpose. Subsequently, we parallelized the code by rewriting the time integrators as graphics card kernels. The GPU accelerated model achieved an amazing 500 times performance boost from the single-processor version. Finally, to make our model comparable to Fortran MPAS-Ocean, we wrote methods to divide the computational work over multiple cores of a cluster with MPI. We then performed equivalent simulations with the Julia and Fortran codes to compare the speeds, and learn how useful Julia might be for climate modeling and other HPC applications.


Currently only includes gravity, coriolis terms (no non-linear terms).



## To replicate, use, or develop:

### Required packages
1. Install the latest version of the Julia language: https://julialang.org/downloads/
2. Install jupyter notebook/lab: https://jupyter.org/install
3. (Optional) if you want to run the distributed (MPI) or graphics card (GPU) simulations on your own computer, install openmpi or CUDA respectively (this will require a MPI-compatible cluster or NVIDIA GPU, respectively). Otherwise, only serial-CPU versions of the simulation can be run.

### Set up
3. Clone the directory and open it <br>
    $ `git clone git@github.com:/robertstrauss/MPAS_Ocean_Julia` or `git clone https://github.com/robertstrauss/MPAS_Ocean_Julia.git` <br>
    $ `cd MPAS_Ocean_Julia` <br>
4. Install required Julia packages (may take some time to download) <br>
    $ `julia` <br>
    julia> `] activate .` <br>
    pkg> `instantiate` <br>
5. (Optional) if you want to run the distributed simulation tests with mpi, install the julia version of the tool for starting a distributed script across nodes: <br>
    julia> `using MPI` <br>
    julia> `MPI.install_mpiexecjl()` <br>
    julia> `exit()` <br>
    a. (Optional) add the tool to your path for easier use: <br>
    $ `ln -s ~/.julia/bin/mpiexecjl <some place on your $PATH>` <br>
    

### Running simulation
There are two main experiments done in this project: <br>
* performance benchmarking of simulation using the graphics card (GPU), compared with the control of simulating on the CPU
* scaling tests run on high-performance clusters, using MPI to distribute simulation across many nodes
As well as verification/validation done to back-up the results.

#### CPU-GPU performance test
Data from performance tests run on Strauss's computer are included in the `./output/asrock/` directory.
To re-run the simulation on your machine and create new data (or just to look at simulation visualization and other information about the simulation):
* (A CUDA-compatible device (computer with an NVIDIA GPU) is required.)
* Start a jupyter server and open the notebook `./GPU_CPU_performance_comparison_meshes.ipynb`.
* Run all the cells to create data files of the performance at `./output/asrock/serialGPU_timing/coastal_kelvinwave/steps_20/resolution_64x64/nvlevels_100/`.

Whether creating new data or using the included, the table from the paper can be recreated:
* Open and run the notebook `./GPU_CPU_performance_comparison_meshes.ipynb` to recreate the tables shown in the paper.

#### MPI (distributed) scaling test
Data from tests run on NERSC's cori-haswell are included in the `./output/` directory.
To re-run the simulations on your cluster and create new data:
* Run `./scaling_test/scaling_test.jl` using mpi, and specify the tests to do <br>
    $ `~/.julia/bin/mpiexecjl --project=. -n <nprocs> julia ./scaling_test/scaling_test.jl <cells> <samples> <proccounts>` <br>
    Where \<nprocs\> = the number of processors needed for the trial distributed over the most processors (should be power of 2), <br>
        \<cells\> = the width of the mesh to use for all trials (128, 256, 512 used in paper), <br>
        \<samples\> = the number of times to re-run a simulation for each trial (reccomended at least 2, the first run is often an outlier), <br>
        \<proccounts\> = a list of how many processor should be used for each trial (use powers of 2, seperated by commas, e.g. 2,8,16,32,256 would run 5 trials, the first distributed over 2 procs, then 8, then 16, etc....) Leave blank to run all powers of 2 up to the number of allocated processors (\<nprocs\>). <br>
* The results will be saved in `./output/kelvinwave/resolution<cells>x<cells>/procs<nprocs>/steps10/nvlevels100/` <br>

The plots from the paper can be recreated in the notebook `./output/kelvinwave/performanceplost.ipynb` using our data or your new data.

#### Verifying/validating the model
To verify the implementation of the mesh, see `./Operator_testing.ipynb`.
To validate the implementation of the shallow water equations, see `KevlinWaveConvergenceTest.ipynb`.
