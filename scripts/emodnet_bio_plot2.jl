using OceanPlot
using NCDatasets
using PyPlot
using Glob
using CSV
using DataFrames
using Proj4
using Dates
using Statistics
using DIVAnd
using DIVAndNN
using PyCall

include("emodnet_bio_grid.jl")

data_analysis = DIVAndNN.Format2020(expanduser("~/tmp/Emodnet-Bio2020/CSV-split"),"analysis")
data_validation = DIVAndNN.Format2020(expanduser("~/tmp/Emodnet-Bio2020/CSV-split"),"validation")


function pyo(a::Array{T,N}) where {T,N}
    numpy_ma = PyCall.pyimport("numpy").ma
    pycall(numpy_ma.array, Any, coalesce.(a,zero(T)), mask=isnan.(a))
end

function plotanalysis(fname)
    ds = NCDataset(fname)
    value_analysis = nomissing(ds["probability"][:,:],NaN)
    sname = split(basename(fname),"_")[2]
    close(ds)

    lon_a,lat_a,obstime_a,value_a,ids_a = DIVAndNN.loadbyname(data_analysis,years,sname)
    lon_cv,lat_cv,obstime_cv,value_cv,ids_cv = DIVAndNN.loadbyname(data_validation,years,sname)

    XY = DIVAnd.ndgrid(gridlon,gridlat)

    mv_a = DIVAndNN.binobs((lon_a,lat_a),value_a,XY)
    mv_cv = DIVAndNN.binobs((lon_cv,lat_cv),value_cv,XY)

    all_data = vcat(value_analysis[:],mv_a[:],mv_cv[:])

    #cmap = PyPlot.cm.hot_r
    cmap = PyPlot.cm.plasma
    orientation = "horizontal"

    function decoration()
        xlim(extrema(gridlon))
        ylim(extrema(gridlat))
        OceanPlot.plotmap()
        OceanPlot.set_aspect_ratio()
    end

    function obsplot(mv,titlestr,cl_prop)
        #=
        sel = value .== 0
        scatter(lon[sel],lat[sel],14,value[sel],cmap = cmap,marker = "x")
        sel = value .== 1
        scatter(lon[sel],lat[sel],4,value[sel],cmap = cmap,marker = "o")
        =#

        pcolor(gridlon,gridlat,pyo(copy(mv')),cmap = cmap)
        clim(cl_prop)
        title(titlestr)
        #colorbar(orientation=orientation)
        decoration()
    end


    fig.suptitle(sname,style="italic")

    subplot(2,2,1);
    imm = pcolor(gridlon,gridlat,value_analysis',cmap = cmap);

    #cl_prop = extrema(value_analysis[isfinite.(value_analysis)]);
    cl_prop = extrema(all_data[isfinite.(all_data)]);

    title("(a) Probability of occurrence");
    clim(cl_prop)
    #colorbar(orientation=orientation)
    decoration()


    cbar_ax = fig.add_axes([0.55, 0.35, 0.35, 0.025]); fig.colorbar(imm, cax=cbar_ax,orientation=orientation)

    subplot(2,2,2);
    #scatter(lon_a,lat_a,10,value_a,cmap = cmap)
    obsplot(mv_a,"(b) Data used in the analysis",cl_prop)

    #@show mean(value_analysis[isfinite.(value_analysis)])
    #@show mean(value_a)
    #@show mean(value_cv)

    subplot(2,2,3);
    obsplot(mv_cv,"(c) Validation data",cl_prop)

    figname = replace(fname,".nc" => ".png")
    savefig(figname)
end

# define
# outdir = ...

fig = figure(figsize = (7,7))
for fname in glob("*nc",outdir)
#for fname in glob("*nc",outdir)[1:1]
    clf();
    #close("all")
    @info(fname)
    plotanalysis(fname)
end
