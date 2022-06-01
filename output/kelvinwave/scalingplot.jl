using DataFrames, CSV

CODE_ROOT = pwd()

include(CODE_ROOT * "/visualization.jl")

nCellsX = 64
nsteps = 10
filename = "2022-05-31T23:39:18.604.txt"
df = DataFrame(CSV.File(CODE_ROOT * "/output/kelvinwave/resolution$(nCellsX)x$(nCellsX)/steps$(nsteps)/$(filename)"))


fortranmeans = [
0.4984766666666667
0.6779433333333333
0.82572
0.6655300000000001
1.1103433333333335
1.4288966666666667
2.1260733333333337
3.2138833333333334
][end:-1:1] ./ 10 # numbers for 10 steps


runs = [:time2, :time3, :time4]
juliameans = 1/length(runs) * sum(Array(df[:,runs]), dims=2) ./ nsteps

fig, ax = plt.subplots()

# ax.loglog(nprocs, fortranmeans, label="fortran")
ax.loglog(df.procs, juliameans, label="julia")
# ax.loglog(df.procs, df.time1, label="compile")
ax.loglog(df.procs, fortranmeans, label="fortran")

ax.set_xticks(df.procs)
ax.set_yticks([0.03, 0.06, 0.12, 0.24])

ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

ax.legend()
ax.set_xlabel("number of processors")
ax.set_ylabel("wallclock time per step (s)")

display(fig)

# fig.savefig("./fortran_x_julia_scaling")
