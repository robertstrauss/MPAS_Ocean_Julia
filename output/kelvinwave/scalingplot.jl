using DataFrames, CSV

CODE_ROOT = pwd()

include(CODE_ROOT * "/visualization.jl")

# nCellsX = 64
# nsteps = 10
# # filename = "2022-05-31T23:39:18.604.txt"
# filename = "2022-06-01T18:00:15.986.txt"
# df = DataFrame(CSV.File(CODE_ROOT * "/output/kelvinwave/resolution$(nCellsX)x$(nCellsX)/steps$(nsteps)/$(filename)"))
# fname = CODE_ROOT * "/output/kelvinwave/resolution64x64/steps10/2022-06-01T18:52:28.303.txt"
nCellsX = 128
fname = CODE_ROOT * "/output/kelvinwave/resolution$(nCellsX)x$(nCellsX)/steps10/2022-06-01T20:17:04.919.txt"
df = DataFrame(CSV.File(fname))

fortranmeans64 = [ 0.0
                   0.0
                  0.4984766666666667
                  0.6779433333333333
                  0.82572
                  0.6655300000000001
                  1.1103433333333335
                  1.4288966666666667
                  2.1260733333333337
                  3.2138833333333334
][end:-1:1] ./ 4 # numbers for 10 steps

fortranmeans128 = [ 0.520183333333333
                  0.68918
                  0.858013333333333
                  1.364753333333333
                  2.378356666666667
                  2.109033333333333
                  3.83726
                  5.057326666666666
                  7.890809999999999
                  13.06165
][end:-1:1]  ./4

fortranmeans = fortranmeans128

runs = [:sim_time2, :sim_time3, :sim_time4, :sim_time5, :sim_time6]
mpis = [:mpi_time3, :mpi_time4, :mpi_time5, :mpi_time6]
juliasimmean = 1/length(runs) * sum(Array(df[:,runs]), dims=2)# ./ nsteps
juliampimean = 1/length(mpis) * sum(Array(df[:,mpis]), dims=2) #./ nsteps
juliameans = juliasimmean .+ juliampimean


perfectjulia = juliameans[1] * df.procs[1] ./ df.procs


fig, ax = plt.subplots()
# ax.loglog(nprocs, fortranmeans, label="fortran")
ax.loglog(df.procs, juliameans, label="julia total")
# for (i, run) in enumerate(runs)
#     ax.loglog(df.procs, Array(df[:,run]), label="julia sim $i")
# end
# for (i, run) in enumerate(mpis)
#     ax.loglog(df.procs, Array(df[:,run]), label="julia mpi $i")
# end
# ax.loglog(df.procs, Array(df[:,runs], label="julia")
# ax.loglog(df.procs, Array(df[:,runs], label="julia")
# ax.loglog(df.procs, Array(df[:,runs], label="julia")
ax.loglog(df.procs, perfectjulia, label="perfect scaling", linestyle=":", color="grey")

ax.loglog(df.procs, juliasimmean, label="julia sim only")
ax.loglog(df.procs, juliampimean, label="julia mpi only")
# ax.loglog(df.procs, df.time1, label="compile")
ax.loglog(df.procs, fortranmeans, label="fortran")


ax.set_xticks(df.procs)
ax.set_yticks([0.03, 0.06, 0.12, 0.24])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xlabel("number of processors")
ax.set_ylabel("wallclock time per step (s)")

ax.set_title("$(nCellsX)x$(nCellsX)")
ax.legend()

display(fig)

fig.savefig("$(fname)_scaling$(nCellsX)x.png")
