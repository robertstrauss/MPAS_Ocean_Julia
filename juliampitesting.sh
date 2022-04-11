#!/bin/sh
#SBATCH --job-name="julia mpas ocean"
#SBATCH --nodes=1
#SBATCH --output=./logs/mpitesting%j.log

scontrol show hostname $SLURM_JOB_NODELIST | perl -ne 'chomb; print "$_"x1'> mpi_hostnames

export PATH="$PATH:/global/homes/r/rstrauss/miconda3/envs/rusty_mpas_julia/bin"

mpiexecjl -n 4 -hostfile mpi_hostnames julia --threads 1 /global/homes/r/rstrauss/repos/MPAS_Ocean_Julia/juliampitesting.jl
