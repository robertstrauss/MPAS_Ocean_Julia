#!/usr/bin/env python
# coding: utf-8

# In[1]:


CODE_ROOT = pwd()


# In[2]:


using DelimitedFiles
include(CODE_ROOT * "/visualization.jl")

function wrap_regex(str::AbstractString, maxlen = 92)
    replace(str, Regex(".{1,$maxlen}( |\$)") => @s_str "\\0\n")
end
# In[11]:


function convergence_plot(operator_name, ncells, ErrorNorm, errortype; decimals=2)
    nCellsx = sqrt.(ncells)
    
    A = [log10.(nCellsx)    ones(length(nCellsx))]
    m, c = A \ log10.(ErrorNorm)
    y = m*log10.(nCellsx) .+ c
    y = 10 .^ y
    
    slopestr ="$(round(m,digits=decimals))"
    while length(split(slopestr, ".")[end]) < decimals
        slopestr *= "0"
    end

    fig, ax = subplots(1,1, figsize=(9,9))
    fig.tight_layout()
    loglog(nCellsx, ErrorNorm, label="$errortype error norm", marker="s", markersize=12.5, linestyle="None", color="black")
    loglog(nCellsx, y, label="Best fit line, slope=$(slopestr)", color="black", linewidth=2)
    ax.set_title("Convergence of $operator_name",
        fontsize=30, fontweight="bold", y=1.035)
    ax.legend(loc="upper right", fontsize=22.5)
    ax.set_xlabel("Number of cells in zonal direction", fontsize=25, labelpad=7.5)
    ax.set_ylabel("$errortype error norm of $(split(lowercase(operator_name), '\n')[1])", fontsize=25, labelpad=7.5)
    ax.tick_params(axis="both", which="both", labelsize=20)

    grid(which="both")
    
    return fig, ax
end


# In[4]:


function readconvergencefile(operator_name, periodicity)
    fname = CODE_ROOT * "/output/operator_convergence/$operator_name/$periodicity/"
    fname = filter(x->x[end-3:end] == ".txt", readdir(fname, join=true))[end]
    return readdlm(fname), fname
end


# In[5]:


operatortitles = Dict{String,String}(
    "Numerical Curl (at Cell Center)" => "Curl Operator \nInterpolated to Cell Centers",
    "Numerical Curl (at Vertex)" => "Curl Operator \nCalculated at Vertices",
    "Numerical Divergence" => "Divergence Operator \nat Cell Centers",
    "Numerical Gradient" => "Gradient Operator \nNormal to Edge",
    "Tangential Velocity" => "Tangential Velocity \nalong Edges"
)


# In[12]:


savefigs = true
for operator_name in keys(operatortitles)
    data, fname = readconvergencefile(operator_name, "Periodic")
    println(fname)
    ncells = data[1,:]
    L2ErrorNorm = data[4,:]
    fig, ax = convergence_plot(operatortitles[operator_name], ncells, L2ErrorNorm, "L\$^2\$")
    display(fig)
    
    if savefigs
        fig.savefig("$(fname)_convergenceplot.pdf", bbox_inches="tight")
        println("saved file at $(fname)_convergenceplot.pdf")
    end
end


# In[ ]:





# In[ ]:




