#!/bin/bash

TESTNAME="17jun-1500"

sbatch -A e3sm -C cpu --ntasks-per-node=64 --qos=regular -o sbatch-${TESTNAME}/16x%A.out  -e sbatch-${TESTNAME}/16x%A/%t.err  -t 2:00  -N 1  mpiexecjl -n 64  julia --project=. ./scaling_test.jl 16  6
sbatch -A e3sm -C cpu --ntasks-per-node=64 --qos=regular -o sbatch-${TESTNAME}/32x%A.out  -e sbatch-${TESTNAME}/32x%A/%t.err  -t 3:00  -N 2  mpiexecjl -n 128 julia --project=. ./scaling_test.jl 32  6
sbatch -A e3sm -C cpu --ntasks-per-node=64 --qos=regular -o sbatch-${TESTNAME}/64x%A.out  -e sbatch-${TESTNAME}/64x%A/%t.err  -t 5:00  -N 8  mpiexecjl -n 512 julia --project=. ./scaling_test.jl 64  6
sbatch -A e3sm -C cpu --ntasks-per-node=64 --qos=regular -o sbatch-${TESTNAME}/128x%A.out -e sbatch-${TESTNAME}/128x%A/%t.err -t 10:00 -N 32 mpiexecjl -n 2048 julia --project=. ./scaling_test.jl 128 6
sbatch -A e3sm -C cpu --ntasks-per-node=64 --qos=regular -o sbatch-${TESTNAME}/256x%A.out -e sbatch-${TESTNAME}/256x%A/%t.err -t 15:00 -N 64 mpiexecjl -n 4096 julia --project=. ./scaling_test.jl 256 6
sbatch -A e3sm -C cpu --ntasks-per-node=64 --qos=regular -o sbatch-${TESTNAME}/512x%A.out -e sbatch-${TESTNAME}/512x%A/%t.err -t 25:00 -N 64 mpiexecjl -n 4096 julia --project=. ./scaling_test.jl 512 6
