using DataFrames, CSV, DelimitedFiles

CODE_ROOT = pwd()

include(CODE_ROOT * "/visualization.jl")

nCellsX = 512
fname = CODE_ROOT * "/output/kelvinwave/resolution$(nCellsX)x$(nCellsX)/procs1024/steps10/nvlevels100/"
fname = filter(x->x[end-3:end] == ".txt", readdir(fname, join=true))[end]
# fname *= "2022-08-24T16:27:20.125.txt"
# fname *= "2022-08-24T16:40:29.640.txt"
# "/output/kelvinwave/resolution$(nCellsX)x$(nCellsX)/steps10/2022-06-01T20:17:04.919.txt"
df = DataFrame(CSV.File(fname))

runs = [:sim_time2, :sim_time3, :sim_time4, :sim_time5, :sim_time6]
mpis = [:mpi_time3, :mpi_time4, :mpi_time5, :mpi_time6]
juliasimmean = 1/length(runs) * sum(Array(df[:,runs]), dims=2)
juliampimean = 1/length(mpis) * sum(Array(df[:,mpis]), dims=2)
juliameans = juliasimmean .+ juliampimean

perfectjulia = juliameans[1] * df.procs[1] ./ df.procs



fortranfname = readdir(CODE_ROOT * "/output/kelvinwave/fortranperformance/resolution$(nCellsX)x$(nCellsX)/", join=true)[1]
fortrantiming = readdlm(fortranfname, skipstart=7)
# julia is using forward-backward, fortran is using RK4, so to make comparison fair we divide fortran times by 4 since its doing about 4x as much work
fortranmeans = fortrantiming[:,end-1] ./ 4
# fortransypd = fortrantiming[:,end] .* 4
fortranprocs = fortrantiming[:,1]
# factor = ( fortransypd .* fortranmeans )[1]



# juliasypd = factor ./ juliameans
# perfectjuliasypd = factor ./ perfectjulia

linewidth = 1
linestyle = "-"
markersize = 10
tickfontsize = 15
labelfontsize = 20

fig, ax = plt.subplots(figsize=(9,9))
ax.loglog(df.procs, juliameans, label="julia", linewidth=linewidth,linestyle="-",color="red",marker="s",markersize=markersize)
ax.loglog(df.procs, perfectjulia, label="perfect scaling", linestyle=":", color="grey")
ax.loglog(fortranprocs, fortranmeans, label="fortran", linewidth=linewidth,linestyle="--",color="blue",marker="D",markersize=markersize)

# ax.loglog(df.procs, Array(df[:,runs], label="julia")
# ax.loglog(df.procs, juliasimmean, label="julia sim only")
# ax.loglog(df.procs, juliampimean, label="julia mpi only")


ax.set_xticks(df.procs)
ax.tick_params(axis="x", labelsize=tickfontsize)
ax.tick_params(axis="y", labelsize=tickfontsize)
# ax.set_yticks([minimum(hcat(juliameans,fortranmeans)),
#         sum(hcat(juliameans,fortranmeans))/(length(juliameans)+length(fortranmeans)),
#         maximum(hcat(juliameans,fortranmeans))])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xlabel("Number of processors", fontsize=labelfontsize)
ax.set_ylabel("Wall clock time elapsed during computation (s)", fontsize=labelfontsize)
ax.set_title("Scaling of Coastal Kelvin Wave Simulation \non Cori-Haswell with $(nCellsX)x$(nCellsX) Hexagonal Mesh", fontsize=21, fontweight="bold")
ax.legend(loc="lower right")

ax.grid(which="both")
plt.tight_layout()

display(fig)

fig.savefig("$(fname)_scaling$(nCellsX)x.png")

fig.savefig("/tmp/plot.png")
