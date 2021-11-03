#!/bin/sh
#SBATCH --job-name="julia mpas ocean"
#SBATCH --nodes=64
#SBATCH --output=./logs/fortrancomparison%j.log


scontrol show hostname $SLURM_JOB_NODELIST | perl -ne 'chomb; print "$_"x1'> mpi_hostnames

export PATH="$PATH:/global/homes/r/rstrauss/miconda3/envs/rusty_mpas_julia/bin"

# jlpath="/global/homes/r/rstrauss/miconda3/envs/Rusty_MPAS_Julia/bin"
# conda init bash
# conda activate Rusty_MPAS_Julia
mpiexecjl -n 4096 -hostfile mpi_hostnames julia --threads 1 /global/homes/r/rstrauss/repos/MPAS_Ocean_Julia/fortran-1000x1000-roughcomparison.jl
