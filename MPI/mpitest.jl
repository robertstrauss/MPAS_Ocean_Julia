#!/global/homes/r/rstrauss/julia-1.6.1/bin/julia
using MPI
MPI.Init()


comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
root = 0
commsize = MPI.Comm_size(comm)

println("rank: $rank, threads: $(Threads.nthreads()), hostname: $(gethostname())")
