#!/bin/sh
#SBATCH --job-name="julia mpas ocean"
#SBATCH --nodes=4
#SBATCH --output=./logs/%j.log


scontrol show hostname $SLURM_JOB_NODELIST | perl -ne 'chomb; print "$_"x1'> mpi_hostnames

jlpath="/global/homes/r/rstrauss/miconda3/envs/rusty_mpas_julia/bin"
$jlpath/mpiexecjl -n 4 -hostfile mpi_hostnames $jlpath/julia ./distributed_simulation_test.jl
