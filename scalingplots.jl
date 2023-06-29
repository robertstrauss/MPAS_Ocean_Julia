#!/usr/bin/env python
# coding: utf-8

# In[4]:


using DataFrames, CSV, DelimitedFiles, Statistics

CODE_ROOT = pwd()

JULIA_DATA_ROOT = CODE_ROOT * "/output"

FORTRAN_DATA_ROOT = CODE_ROOT * "/data/new-fortran-timing/"

include(CODE_ROOT * "/visualization.jl")


# In[5]:


# the maximum number of processes each mesh resolution was parallelized with
# each process should have at least 4 cells (nCellsX^2)/maxprocs[nCellsX] >= 4
# julia and fortran run with different limits due to slightly different mesh partitionings
maxprocs = Dict{Int64, Int64}(
    16 => 64,
    32 => 128,
    64 => 512,
    128 => 2048,
    256 => 4096,
    512 => 4096,
    )

ftprocs = Dict{Int64, Vector{Tuple{Int64, Int64}}}(
    16 => [(1, 64)],
    32 => [(1,64)],
    64 => [(1,64)],
    128 => [(64, 4096), (1, 256)],
    256 => [(64, 4096), (1, 256)],
    512 => [(64, 4096), (1, 256)],
    )
("")


# In[6]:


function latestfile(dir, filterfunc)
    fnames = filter(filterfunc, readdir(dir, join=true))
    @assert length(fnames) > 0 "no file in $dir matching pattern!"
    return fnames[end]
end


# In[19]:


function juliatimes(nCellsX, jprocs=4096)
    # read julia data from file
    fname = latestfile(JULIA_DATA_ROOT * "/kelvinwave/resolution$(nCellsX)x$(nCellsX)/procs$(jprocs)/steps10/nvlevels100/", x->x[end-3:end] == ".txt")
    df = DataFrame(CSV.File(fname))
    
    # get the simulation time and communication time columns, excluding the first one (this is usually slower from warm-up)
    simkeys = filter(col->startswith(col,"sim_time"), names(df))[2:end]
    mpikeys = filter(col->startswith(col,"mpi_time"), names(df))[2:end]
    
    simtimes = Array(df[:,simkeys])
    mpitimes = Array(df[:,mpikeys])
    runtimes = simtimes + mpitimes
    
    juliarunmeans = dropdims( mean(runtimes, dims=2), dims=2)
    juliampimeans = dropdims( mean(mpitimes, dims=2), dims=2)
    
    juliarunminmax = permutedims(hcat(minimum(runtimes, dims=2), maximum(runtimes, dims=2)), (2,1))
    juliampiminmax = permutedims(hcat(minimum(mpitimes, dims=2), maximum(mpitimes, dims=2)), (2,1))
    
    return juliarunmeans, juliampimeans, df.procs, juliarunminmax, juliampiminmax, fname
end
function fortrantimes(nCellsX, fprocs=[(64, 4096), (1, 256)], plane=64)
    fortrandir = FORTRAN_DATA_ROOT * "resolution$(nCellsX)x$(nCellsX)"
    
    # read fortran data from files
    runtimes = Dict{Int64, Float64}()
    mpitimes = Dict{Int64, Float64}()
    runminmax = Dict{Int64, Tuple{Float64, Float64}}()
    mpiminmax = Dict{Int64, Tuple{Float64, Float64}}()
    for (minprocs, maxprocs) in fprocs
        fortranfnamerun = latestfile(fortrandir,
            x -> occursin("$(maxprocs)pmax", x) && occursin("$(minprocs)pmin", x) && occursin("$(plane)plane", x) && occursin("runtime", x))
        fortranfnamempi = latestfile(fortrandir,
            x -> occursin("$(maxprocs)pmax", x) && occursin("$(minprocs)pmin", x) && occursin("$(plane)plane", x) && occursin("halotime", x))
        fortranruntiming = readdlm(fortranfnamerun, skipstart=8)
        fortranmpitiming = readdlm(fortranfnamempi, skipstart=8)
        for (i, p) in enumerate(fortranruntiming[:,1])
            runtimes[p] = fortranruntiming[i,end]
            mpitimes[p] = fortranmpitiming[i,end]
            runminmax[p] = (minimum(fortranruntiming[i,2:end]), maximum(fortranruntiming[i,2:end]))
            mpiminmax[p] = (minimum(fortranmpitiming[i,2:end]), maximum(fortranmpitiming[i,2:end]))
        end
    end
    procs = sort(collect(keys(runtimes)))
    fruntimes = collect([runtimes[p] for p in procs])
    fmpitimes = collect([mpitimes[p] for p in procs])
    frunminmax = collect([runminmax[p][k] for k in 1:2, p in procs])
    fmpiminmax = collect([mpiminmax[p][k] for k in 1:2, p in procs])
    return fruntimes, fmpitimes, procs, frunminmax, fmpiminmax, fortrandir
end
function juliafortrantimesplits(nCellsX, jprocs="max", fprocs="max")
    if jprocs == "max"
        jprocs = maxprocs[nCellsX]
    end
    if fprocs == "max"
        fprocs = ftprocs[nCellsX]
    end
    juliarunmean, juliampimean, juliaprocs, juliarunminmax, juliampiminmax, juliafname = juliatimes(nCellsX, jprocs)
    fortranrunmean, fortranmpimean, fortranprocs, fortranrunminmax, fortranmpiminmax, fortranfnamerun = fortrantimes(nCellsX, fprocs)
    return juliarunmean, juliampimean, juliaprocs, fortranrunmean, fortranmpimean, fortranprocs, juliafname, fortranfnamerun, juliarunminmax, juliampiminmax, fortranrunminmax, fortranmpiminmax
end


# In[18]:

# 
# for nCellsX in [16, 32, 64, 128, 256, 512]
#     frun, fmpi, fprocs, frunminmax, fmpiminmax, ffile = fortrantimes(nCellsX, ftprocs[nCellsX], 64)
#     
#     println("$(nCellsX)x, 64-plane: any mpi times > total? $(any(fmpi .> frun))")
#     if any(fmpi .> frun)
#         println("file: ", ffile)
#         println("total time: ", frun)
#         println("mpi time: ", fmpi)
#     else
#         println("minimum computation time (total-mpi): ", minimum(frun - fmpi))
#     end
# end


# In[20]:


linewidth = 1
linestyle = "-"
markersize = 10
tickfontsize = 15
labelfontsize = 22.5
titlefontsize = 25
blue = "blue"
red = "red"


# In[21]:


function strongscalingplot(juliameans, juliaprocs, fortranmeans, fortranprocs, nCellsX, juliaerr=0, fortranerr=0; perfect=true)

    perfectjulia = juliameans[1] * juliaprocs[1] ./ juliaprocs
    
    fig, ax = plt.subplots(figsize=(9,9))
    ax.set_xscale("log")
    ax.set_yscale("log")
    
    ax.errorbar(juliaprocs, juliameans, yerr=juliaerr, label="Julia", linewidth=linewidth,linestyle="-",
                capsize=5, color=red,marker="s",markersize=markersize)
    ax.errorbar(fortranprocs, fortranmeans, yerr=fortranerr, label="Fortran", linewidth=linewidth,
                capsize=5, linestyle="--", color=blue,marker="D",markersize=markersize)
    if perfect
        ax.plot(juliaprocs, perfectjulia, label="Perfect scaling", linestyle=":", color="black", linewidth=2)
    end

    ax.set_xticks(juliaprocs)
    ax.tick_params(axis="x", labelsize=tickfontsize)
    ax.tick_params(axis="y", labelsize=tickfontsize)
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    # ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    
    ax.set_xlabel("Number of processors", fontsize=labelfontsize, labelpad=10)
    ax.set_ylabel("Wall-clock time (s)", fontsize=labelfontsize, labelpad=10)
    ax.set_title("$(nCellsX)x$(nCellsX) Hexagonal Mesh", fontsize=titlefontsize, fontweight="bold", y=1.02)
    ax.legend(loc="upper right", fontsize=labelfontsize-2.5)

    ax.grid(which="both")
    plt.tight_layout()

    return fig, ax
end


# In[25]:


savefigs=true


# In[23]:

# 
# # for plane in [32, 64]
# for nCellsX in [128, 256, 512]
#     jrun, jmpi, jprocs, jrunminmax, jmpiminmax, jfile = juliatimes(nCellsX, maxprocs[nCellsX])
#     frun32, fmpi32, fprocs32, frunminmax32, fmpiminmax32, _ = fortrantimes(nCellsX, [(64, 4096)], 32)
#     frun64, fmpi64, fprocs64, frunminmax64, fmpiminmax64, _ = fortrantimes(nCellsX, ftprocs[nCellsX], 64)
# 
#     fig, ax = strongscalingplot(jrun, jprocs, frun64, fprocs64, nCellsX,
#                         abs.(jrun' .- jrunminmax), abs.(frun64' .- frunminmax64))
#     ax.errorbar(fprocs32, frun32, yerr=abs.(frun32' .- frunminmax32), label="Fortran plane 32", linewidth=linewidth,
#                 capsize=5, linestyle="--", color="cyan", marker="D",markersize=markersize)
#     ax.legend()
#     # ax.set_title("Plane=$(plane) $(nCellsX)x$(nCellsX) Hexagonal Mesh\nComputation and Communication", fontsize=titlefontsize, fontweight="bold", y=1.02)
# end
# # end


# In[26]:


for nCellsX in [128, 256, 512] # 16, 32, 64, 
    jrun, jmpi, jprocs, jrunminmax, jmpiminmax, jfile = juliatimes(nCellsX, maxprocs[nCellsX])
    frun, fmpi, fprocs, frunminmax, fmpiminmax, ffile = fortrantimes(nCellsX, ftprocs[nCellsX])
    
    # only take the data points they both have
    jmatch = findall(i -> i in fprocs, jprocs)
    fmatch = findall(i -> i in jprocs, fprocs)
    frun, fmpi, fprocs, frunminmax, fmpiminmax = frun[fmatch], fmpi[fmatch], fprocs[fmatch], frunminmax[:,fmatch], fmpiminmax[:,fmatch]
    jrun, jmpi, jprocs, jrunminmax, jmpiminmax = jrun[jmatch], jmpi[jmatch], jprocs[jmatch], jrunminmax[:,jmatch], jmpiminmax[:,jmatch]
    
    fig, ax = strongscalingplot(jrun, jprocs, frun, fprocs, nCellsX,
                        abs.(jrun' .- jrunminmax), abs.(frun' .- frunminmax))
    ax.set_title("$(nCellsX)x$(nCellsX) Hexagonal Mesh\nComputation and Communication", fontsize=titlefontsize, fontweight="bold", y=1.02)
    if savefigs
        savename = replace(jfile[length(CODE_ROOT)+2:end], ("/" => "--"))
        fig.savefig("$(CODE_ROOT)/plots/strong_scaling/total_$(savename).pdf", bbox_inches="tight")
        println("saved at $(CODE_ROOT)/plots/strong_scaling/total_$(savename).pdf")
    end
    
    ## gpu
    # gpufname = latestfile(CODE_ROOT *  "/output/inertiagravitywave/NVIDIA A100-SXM4-40GB/resolution$(nCellsX)x$(nCellsX)/steps10/nvlevels100/", x -> x[end-3:end] == ".txt")
    # gputimes = readdlm(gpufname)
    # simtimes = gputimes[1,:]
    # bustimes = gputimes[2,:]
    # meangpu = mean(simtimes[2:end])
    # meangpucom = mean(bustimes[2:end])
    # plt.axhline(y=meangpu, color="g")
    # plt.axhline(y=meangpucom, color="g", linestyle="--")

    
    fig1, ax1 = strongscalingplot(jrun - jmpi, jprocs, frun - fmpi, fprocs, nCellsX,
                            abs.((jrun' - jmpi') .- (jrunminmax .- jmpiminmax)), abs.((frun' - fmpi') .- (frunminmax .- fmpiminmax)))
    ax1.set_title("$(nCellsX)x$(nCellsX) Hexagonal Mesh\nComputation Only", fontsize=titlefontsize, fontweight="bold", y=1.02)
    
    # exclude communication time at 1 process, this is meaningless
    jcut1 = findall(p -> p != 1, jprocs)
    fcut1 = findall(p -> p != 1, fprocs)
    fig2, ax2 = strongscalingplot(jmpi[jcut1], jprocs[jcut1], fmpi[fcut1], fprocs[fcut1], nCellsX, 
                            abs.(jmpi' .- jmpiminmax)[:,jcut1], abs.(fmpi' .- fmpiminmax)[:,fcut1], perfect=false)
    ax2.set_title("$(nCellsX)x$(nCellsX) Hexagonal Mesh\nCommunication Only", fontsize=titlefontsize, fontweight="bold", y=1.02)
    
    # fig, ax = strongscalingplot2(jrun .- jmpi, jrunminmax, jmpi, jmpiminmax, jprocs, frun .- fmpi, frunminmax, fmpi, fmpiminmax, fprocs, nCellsX)
    
    if savefigs
        savename = replace(jfile[length(CODE_ROOT)+2:end], ("/" => "--"))
        fig1.savefig("$(CODE_ROOT)/plots/strong_scaling/calc_only_$(savename).pdf", bbox_inches="tight")
        fig2.savefig("$(CODE_ROOT)/plots/strong_scaling/comm_only_$(savename).pdf", bbox_inches="tight")
        println("saved at $(CODE_ROOT)/plots/strong_scaling/calc_only_$(savename).pdf")
        println("saved at $(CODE_ROOT)/plots/strong_scaling/comm_only_$(savename).pdf")
    end
    
    println(jfile[end-60:end])
    println(ffile[end-60:end])
end

function strongscalingplot2(juliasim, juliasimminmax, juliampi, juliampiminmax, juliaprocs, fortransim, fortransimminmax, fortranmpi, fortranmpiminmax, fortranprocs, nCellsX)

    perfectjulia = (juliampi[1] + juliasim[1]) * juliaprocs[1] ./ juliaprocs

    linewidth = 1
    linestyle = "-"
    markersize = 10
    tickfontsize = 15
    labelfontsize = 22.5
    titlefontsize = 25
    blue = "blue"
    red = "red"

    fig, ax = plt.subplots(figsize=(9,9))
    ax.set_xscale("log")
    ax.set_yscale("log")
    
    ax.errorbar(juliaprocs, juliampi, yerr=abs.(juliampi' .-juliampiminmax), label="Julia communication", linewidth=linewidth,linestyle="--",color=red,marker="P",markersize=markersize, capsize=5)
    ax.errorbar(fortranprocs, fortranmpi, yerr=abs.(fortranmpi' .-fortranmpiminmax), label="Fortran communication", linewidth=linewidth,linestyle="--",color=blue,marker="X",markersize=markersize, capsize=5)
    
    ax.errorbar(juliaprocs, juliasim, yerr=abs.(juliasim' .-juliasimminmax), label="Julia computation", linewidth=linewidth,linestyle="-",color=red,marker="s",markersize=markersize, capsize=5)
    ax.errorbar(fortranprocs, fortransim, yerr=abs.(fortransim' .-fortransimminmax), label="Fortran computation", linewidth=linewidth,linestyle="-",color=blue,marker="D",markersize=markersize, capsize=5)
    
    
    ax.loglog(juliaprocs, perfectjulia, label="Perfect scaling", linestyle=":", color="black", linewidth=2)


    ax.set_xticks(juliaprocs)
    ax.tick_params(axis="x", labelsize=tickfontsize)
    ax.tick_params(axis="y", labelsize=tickfontsize)
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    
    ax.set_xlabel("Number of processors", fontsize=labelfontsize, labelpad=10)
    ax.set_ylabel("Wall-clock time (s)", fontsize=labelfontsize, labelpad=10)
    ax.set_title("$(nCellsX)x$(nCellsX) Hexagonal Mesh", fontsize=titlefontsize, fontweight="bold", y=1.02)
    ax.legend(loc="upper right", fontsize=labelfontsize-2.5)

    ax.grid(which="both")
    plt.tight_layout()

    return fig, ax
end
# In[109]:


function timesplitplot(nprocs, comptime, mpitime, nCellsX, info="")
    fontscale = 2
    
    fig, ax = plt.subplots(1,1, figsize=(length(comptime)*3,8))
#     df = DataFrame(processors=nprocs, computation=comptime, communication=mpitime)
#     df.plot(kind="bar", 
#             stacked=true, 
#             colormap="tab10", 
#             figsize=(10, 6))
    totals = comptime + mpitime
    
    ax.bar(string.(Int.(nprocs)), comptime ./ totals, color="blue", label="Computation")
    ax.bar(string.(Int.(nprocs)), mpitime ./ totals, bottom=comptime ./ totals, color="red", label="Communication")
    
    ax.legend(loc="center left", fontsize=20*fontscale, bbox_to_anchor=(1,0.5))
    
    ax.set_yticks([0.0, 0.5, 1.0])
    ax.set_xlabel("Number of processors", fontsize=25*fontscale, labelpad=10.0*fontscale)
    ax.set_ylabel("Proportion of time", fontsize=25*fontscale, labelpad=10.0*fontscale)
    ax.tick_params(axis="x", labelsize=20*fontscale)
    ax.tick_params(axis="y", labelsize=20*fontscale)
    
    # " using $(nCellsX)x$(nCellsX) Hexagonal Mesh"
    ax.set_title("$info: Proportion of Simulation Time \nSpent on Computation and Communication", fontweight="bold", fontsize=30*fontscale, y=1.08)
    
#     plt.tick_params(top="off", bottom="off", left="off", right="off", labelleft="off", labelbottom="on")
    
    # plt.tight_layout()
    
    return fig, ax
end


# In[111]:


nCellsX = 512
juliasim, juliampi, juliaprocs, fortransim, fortranmpi, fortranprocs, juliafname, fortranfname, juliasimminmax, juliampiminmax = juliafortrantimesplits(nCellsX)

fig1, ax = timesplitplot(juliaprocs, juliasim, juliampi, nCellsX, "Julia")
fig2, ax = timesplitplot(fortranprocs[end:-1:1], fortransim[end:-1:1], fortranmpi[end:-1:1], nCellsX, "Fortran")

if savefigs
    savenamejl = replace(juliafname[length(CODE_ROOT)+2:end], ("/" => "--"))
    savenameft = replace(fortranfname[length(CODE_ROOT)+2:end], ("/" => "--"))
    fig1.savefig("$(CODE_ROOT)/plots/time_proportion/julia_$(savenamejl).pdf", bbox_inches="tight")
    # fig1.savefig("$(juliafname)_proportion_sim_mpi.pdf", bbox_inches="tight")
    fig2.savefig("$(CODE_ROOT)/plots/time_proportion/fortran_$(savenameft).pdf", bbox_inches="tight")
end


# In[8]:


function weakscalingplot(;which="total", resolutions = [16, 32, 64, 128, 256, 512], constlines = [64], mode="weak_scaling", gpu=false)
    function constline(ncells, procs)
        if mode == "weak_scaling"
            return ncells ./ procs
        elseif mode == "constant_nprocs"
            return procs
        else 
            error("unknown mode")
        end
    end
    
    # if mode == "constant nprocs" && constlines == [64]
    #     global FORTRAN_DATA_ROOT = CODE_ROOT * "/data/fortran-timing/timing64_plane"
    # else
    #     global FORTRAN_DATA_ROOT = CODE_ROOT * "/data/fortran-timing"
    # end
    
    jltimes   = zeros((  length(resolutions),length(constlines)))
    jlminmaxs = zeros((2,length(resolutions),length(constlines)))
    fttimes   = zeros((  length(resolutions),length(constlines)))
    ftminmaxs = zeros((2,length(resolutions),length(constlines)))
    gputimes  = zeros((  length(resolutions),length(constlines)))
    gpuminmaxs= zeros((2,length(resolutions),length(constlines)))
    fnames = Vector{String}(undef, length(resolutions))
    
    for (i, nCellsX) in enumerate(resolutions)
    
        jrun, jmpi, jprocs, jrunminmax, jmpiminmax, fnames[i] = juliatimes(nCellsX, maxprocs[nCellsX])
        ind = []
        for con in constlines
            ind = vcat(ind, findall(x->x==con, constline(nCellsX^2, jprocs)))
        end
        
        if which == "total"
            jltimes[i,:] = jrun[ind]
            jlminmaxs[:,i,:] = jrunminmax[:,ind]
        elseif which == "mpi"
            jltimes[i,:] = jmpi[ind]
            jlminmaxs[:,i,:] = jrunminmax[:,ind]
        elseif which == "comp"
            jltimes[i,:] = jrun[ind] - jmpi[ind]
            jlminmaxs[:,i,:] = jrunminmax[:,ind] - jmpiminmax[:,ind]
        end
        
        if mode == "constant nprocs" && constlines == [64]
            if nCellsX == 512
                frun, fmpi, fprocs, frunminmax, _, _ = fortrantimes(nCellsX, [(64, 4096)], 64)
            else
                frun, fmpi, fprocs, frunminmax, _, _ = fortrantimes(nCellsX, [(1,256)], 64)
            end
        else
            frun, fmpi, fprocs, frunminmax, fmpiminmax, _ = fortrantimes(nCellsX, ftprocs[nCellsX])
        end
        ind = []
        for con in constlines
            ind = vcat(ind, findall(x->x==con, constline(nCellsX^2, fprocs)))
        end
        
        if which == "total"
            fttimes[i,:] = frun[ind]
            ftminmaxs[:,i,:]  = frunminmax[:,ind]
        elseif which == "mpi"
            fttimes[i,:] = fmpi[ind]
            ftminmaxs[:,i,:]  = fmpiminmax[:,ind]
        elseif which == "comp"
            fttimes[i,:] = frun[ind] - fmpi[ind]
            ftminmaxs[:,i,:] = frunminmax[:,ind] - fmpiminmax[:,ind]
        end
        
        if gpu
            gpufname = latestfile(CODE_ROOT *  "/output/inertiagravitywave/NVIDIA A100-SXM4-40GB/resolution$(nCellsX)x$(nCellsX)/steps10/nvlevels100/", x -> x[end-3:end] == ".txt")
            gputtimes = readdlm(gpufname)
            simtimes = gputtimes[1,:]
            bustimes = gputtimes[2,:]
            if which == "total"
                times = simtimes[2:end] + bustimes[2:end]
            elseif which == "mpi"
                times = bustimes[2:end]
            elseif which == "comp"
                times = simtimes[2:end]
            end
            gpuminmaxs[:,i,:] = vcat(minimum(times), maximum(times))
            gputimes[i,:] .= mean(times)
        end
        
    end
    
    # plot
    println(fnames)
    for (i, con) in enumerate(constlines)
        fig, ax = plt.subplots(1,1, figsize=(9,9))
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylabel("Wall-clock time (s)", fontsize=labelfontsize, labelpad=10)
        
        if which == "total"
            whichname = "Computation and Communication"
        elseif which == "mpi"
            whichname = "Communication Only"
        elseif which == "comp"
            whichname = "Computation Only"
        end
        if mode == "weak_scaling"
            ax.set_title("Weak Scaling: $(con) Cells/Process\n$(whichname)", fontsize=titlefontsize, fontweight="bold", y=1.02)
            ax.set_xlabel("Number of processors", fontsize=labelfontsize, labelpad=10)
            xaxis = Int.(resolutions .^ 2 / con) # = num proccessors
            xlabels = xaxis
        elseif mode == "constant_nprocs"
            ax.set_title("$(con) MPI Processes versus 1 GPU\n$(whichname)", fontsize=titlefontsize, fontweight="bold", y=1.02)
            ax.set_xlabel("Mesh size", fontsize=labelfontsize, labelpad=10)
            xaxis = resolutions .^ 2 # mesh size
            xlabels = string.(resolutions) .* "x" .* string.(resolutions)
        end
        
        julialines = ax.errorbar(xaxis, jltimes[:,i], yerr=abs.(jltimes[:,i]' .-jlminmaxs[:,:,i]), label="Julia MPI", 
                        capsize=5, linewidth=linewidth, linestyle="-", marker="s", markersize=markersize, color=red)
        fortranlines = ax.errorbar(xaxis, fttimes[:,i], yerr=abs.(fttimes[:,i]' .-ftminmaxs[:,:,i]), label="Fortran MPI",
                        capsize=5, linewidth=linewidth, linestyle="--", marker="D", markersize=markersize, color=blue)
        if gpu
            gpuline = ax.errorbar(xaxis, gputimes[:,i], yerr=abs.(gputimes[:,i]' .- gpuminmaxs[:,:,i]), label="Julia GPU",
                        capsize=5, linewidth=linewidth, linestyle=":", marker="o", markersize=markersize, color="green")
        end
        
        low = minimum(vcat(jltimes[:,i], fttimes[:,i]))
        high = maximum(vcat(jltimes[:,i], fttimes[:,i]))
        if gpu
            low = minimum(vcat(low, gputimes))
            high = maximum(vcat(high, gputimes))
        end
        ax.set_ylim(low/4+1e-9, high*4)
        ax.set_xticks(xaxis)
        ax.tick_params(axis="x", labelsize=tickfontsize)
        ax.tick_params(axis="y", labelsize=tickfontsize)
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.set_xticklabels(xlabels)
        # ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.legend(loc="upper center", fontsize=labelfontsize-2.5)
        ax.grid(which="both")
        plt.tight_layout()
        
        if mode == "weak_scaling" && which == "mpi" && xaxis[1] == 1
            # clip out single-processor communication time as this is meaningless
            julialines[1].set_data(xaxis[2:end], jltimes[2:end,i])
            fortranlines[1].set_data(xaxis[2:end], fttimes[2:end,i])
            for errorbar in [ julialines[3][1], fortranlines[3][1] ]
                errorbar.set_alpha(vcat(0, [1 for i in 2:length(xaxis)]))
            end
            for cap in (julialines[2]..., fortranlines[2]...)
                capdata = cap.get_data()
                cap.set_data((capdata[1][2:end], capdata[2][2:end]))
            end
            low = minimum(vcat(jltimes[2:end,i], fttimes[2:end,i]))
            high = maximum(vcat(jltimes[2:end,i], fttimes[2:end,i]))
            if gpu
                low = minimum(vcat(low, gputimes))
                high = maximum(vcat(high, gputimes))
            end
            ax.set_ylim(low/4+1e-9, high*4)
        end
        
        
        if savefigs
            savename = replace(fnames[end][length(CODE_ROOT)+2:end], ("/" => "--"))
            hasgpu = gpu ? "and-gpu_" : ""
            savepath = "$(CODE_ROOT)/plots/$(mode)/$(which)_$(con)_$(hasgpu)$(resolutions[1])x-$(resolutions[end])x_$(savename).pdf"
            fig.savefig(savepath, bbox_inches="tight")
            println("saved at $(savepath)")
        end
    end
end


# In[9]:


# savefigs=false


# In[10]:


for which in ["total", "mpi", "comp"]
    weakscalingplot(which = which, resolutions = [16, 32, 64, 128, 256, 512], constlines = [64, 128, 256])
end


# In[ ]:





# In[64]:


for which in ["total", "mpi", "comp"]
    weakscalingplot(which=which, resolutions = [16, 32, 64, 128, 256, 512], constlines = [64,], mode="constant_nprocs", gpu=true)
end

#=
# In[39]:


fortrantimes(32, ftprocs[32])


# In[ ]:





# In[229]:


nCellsX = 256
jprocs = 4096
fprocs = 4096
fname = latestfile(JULIA_DATA_ROOT * "/kelvinwave/resolution$(nCellsX)x$(nCellsX)/procs$(jprocs)/steps10/nvlevels100/", x->x[end-3:end] == ".txt")
df = DataFrame(CSV.File(fname))

runs = filter(col->startswith(col,"sim_time"), names(df))[2:end]
mpis = filter(col->startswith(col,"mpi_time"), names(df))[2:end]
juliasimmean = 1/length(runs) * sum(Array(df[:,runs]), dims=2)
juliampimean = 1/length(mpis) * sum(Array(df[:,mpis]), dims=2)
juliameans = juliasimmean .+ juliampimean


fortranfnamesim = latestfile(FORTRAN_DATA_ROOT, x -> occursin("$fprocs", x) && occursin("$(nCellsX)x$(nCellsX)", x) && x[end-3:end] == ".txt" && occursin("runtime", x))
fortranfnamempi = latestfile(FORTRAN_DATA_ROOT, x -> occursin("$fprocs", x) && occursin("$(nCellsX)x$(nCellsX)", x)  && x[end-3:end] == ".txt" && occursin("halotime", x))
fortransimtiming = readdlm(fortranfnamesim, skipstart=8)
fortranmpitiming = readdlm(fortranfnamempi, skipstart=8)
fortransimmean = fortransimtiming[:,end] #./ 4
fortranmpimean = fortranmpitiming[:,end]
fortransimstd = dropdims( std(Array(fortransimtiming[:,2:end-1]), dims=2), dims=2)
fortranmpistd = dropdims( std(Array(fortranmpitiming[:,2:end-1]), dims=2), dims=2)
fortranprocs = fortransimtiming[:,1]

juliasimstd = std(Array(df[:,runs]), dims=2)
juliampistd = std(Array(df[:,mpis]), dims=2)

df


# In[59]:


DataFrame(fortransimtiming[end:-1:1,:], :auto)


# In[ ]:


fig, ax = plt.subplots(1,1, figsize=(9,9))

ax.set_xscale("log")
ax.set_yscale("log")

for i in 1:length(mpis)
    mpi = mpis[i]
    ax.loglog(df.procs *0.95, df[:,mpi],               marker="s", linestyle="none", color="red",  markersize=2)
    ax.loglog(fortranmpitiming[:,1], fortranmpitiming[:,i+1], marker="D", linestyle="none", color="blue", markersize=2)
end

ax.loglog(df.procs, juliampimean, marker="D", color="orange")
ax.loglog(fortranmpitiming[:,1], fortranmpimean, marker="D", color="cyan")

ax.set_xlabel("number of processors")
ax.set_ylabel("commmunication time")
("")


# In[ ]:





# In[227]:


nCellsX = 64
gpufname = latestfile(CODE_ROOT *  "/output/inertiagravitywave/NVIDIA A100-SXM4-40GB/resolution$(nCellsX)x$(nCellsX)/steps10/nvlevels100/", x -> x[end-3:end] == ".txt")
gputtimes = readdlm(gpufname)


# In[235]:


jr, jc, jp, _ = juliatimes(64, 512)


# In[236]:


jr


# In[237]:


jp


# In[238]:


jc


# In[ ]:


=#
