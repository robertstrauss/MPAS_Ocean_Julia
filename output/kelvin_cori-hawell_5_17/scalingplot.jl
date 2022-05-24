fig, ax = plt.subplots()

fortrameans = [
 0.82572
 0.66553
 1.11034
 1.42889
 2.12607
 3.21388] ./ 10 # numbers for 10 steps

juliameans = [
 0.7826035610000001
 1.147424507
 1.646721549
 2.509213951
 3.9189539279999996
 6.997029351] ./ 25 # numbers for 25 steps

nprocs = [
 64
 32
 16
  8
  4
  2]


ax.loglog(nprocs, fortranmeans, label="fortran")
ax.loglog(nprocs, juliameans, label="julia")

ax.set_xticks(nprocs)
ax.set_yticks([0.03, 0.06, 0.12, 0.24])

ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

ax.legend()
ax.set_xlabel("number of processors")
ax.set_ylabel("wallclock time per step (s)")

display(fig)

fig.savefig("./fortran_x_julia_scaling")
