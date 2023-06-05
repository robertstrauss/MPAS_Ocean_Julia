#!/bin/bash

sbatch -A e3sm -C cpu --distribution=plane=32 --qos=regular -o sbatch-jun-5-%A.out -e sbatch-jun-5-%A.err -t 2:00 -n 256 mpiexecjl -n 256 julia --project=. ./scaling_test/scaling_test.jl 16 6
sbatch -A e3sm -C cpu --distribution=plane=32 --qos=regular -o sbatch-jun-5-%A.out -e sbatch-jun-5-%A.err -t 2:00 -n 2048 mpiexecjl -n 2048 julia --project=. ./scaling_test/scaling_test.jl 32 6
sbatch -A e3sm -C cpu --distribution=plane=32 --qos=regular -o sbatch-jun-5-%A.out -e sbatch-jun-5-%A.err -t 3:00 -n 4096 mpiexecjl -n 4096 julia --project=. ./scaling_test/scaling_test.jl 64 6
sbatch -A e3sm -C cpu --distribution=plane=32 --qos=regular -o sbatch-jun-5-%A.out -e sbatch-jun-5-%A.err -t 5:00 -n 4096 mpiexecjl -n 4096 julia --project=. ./scaling_test/scaling_test.jl 128 6
sbatch -A e3sm -C cpu --distribution=plane=32 --qos=regular -o sbatch-jun-5-%A.out -e sbatch-jun-5-%A.err -t 10:00 -n 4096 mpiexecjl -n 4096 julia --project=. ./scaling_test/scaling_test.jl 256 6
sbatch -A e3sm -C cpu --distribution=plane=32 --qos=regular -o sbatch-jun-5-%A.out -e sbatch-jun-5-%A.err -t 15:00 -n 4096 mpiexecjl -n 4096 julia --project=. ./scaling_test/scaling_test.jl 512 6
