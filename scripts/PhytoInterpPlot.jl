using PyCall
using PyPlot
using Proj4
using DIVAnd

ccrs = pyimport("cartopy.crs")
gridliner = pyimport("cartopy.mpl.gridliner")
cfeature = pyimport("cartopy.feature")
mticker = pyimport("matplotlib.ticker")
myproj = ccrs.PlateCarree()
coast = cfeature.GSHHSFeature(scale="full");
mpl = pyimport("matplotlib");
cartopyticker = pyimport("cartopy.mpl.ticker")
lon_formatter = cartopyticker.LongitudeFormatter()
lat_formatter = cartopyticker.LatitudeFormatter()

"""
    make_plot_presence_absence(lon, lat, figname)

Plot the locations of absence and presence in 2 subplots and save it as `figname`.

## Examples
```julia-repl
julia> make_plot_presence_absence(lon, lat, "data_locations.png")
```
"""
function make_plot_presence_absence(lon::Vector, lat::Vector, figname::String="")
    data_presence = occurs .== 1;
    data_absence = .!(data_presence);

    fig = PyPlot.figure(figsize=(8, 8))
    ax = PyPlot.subplot(211, projection=myproj)
    ax.plot(lon[data_presence], lat[data_presence], "ko", markersize=.2)
    title("Presence")
    decorate_map(ax)

    ax = PyPlot.subplot(212, projection=myproj)
    ax.plot(lon[data_absence], lat[data_absence], "ko", markersize=.2)
    title("Absence")
    decorate_map(ax)
    if length(figname) > 0
        PyPlot.savefig(figname, dpi=300, bbox_inches="tight")
    end

    PyPlot.close()
end

"""
    decorate_map(ax)

Decorate the map by adding ticks and coastline.

## Examples
```julia-repl
julia> decorate_map(ax)
```
"""
function decorate_map(ax)
    PyPlot.grid(linewidth=0.2)
    ax.add_feature(coast, color=".6",
        edgecolor="k", zorder=5)

    ax.set_xlim(-2., 9.)
    ax.set_ylim(51., 56)
    ax.set_xticks(-2:2.:9.)
    ax.set_yticks(51.:1:56.)
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

end

"""
    plot_heatmap(longrid, latgrid, dens, lonobs, latobs, occurs, titletext, figname)

Plot the heatmap `dens` (computed by `DIVAnd`) on the grid defined by `longrid`, `latgrid`
and add the locations of the observations defined by `lonobs`, `latobs` and `occurs`.

The vector `occurs` (0's and 1's) is used to separate 'presence' and 'absence' data points.

## Examples
```julia-repl
julia> plot_heatmap(longrid, latgrid, dens, lonobs, latobs, occurs,
    "Actinocyclus", "Actinocyclus_data.png")
```
"""
function plot_heatmap(longrid::StepRangeLen, latgrid::StepRangeLen,
    dens::Array, lonobs::Vector, latobs::Vector, occurs::Vector,
    titletext::String, figname::String="")

    data_presence = occurs .== 1;
    data_absence = occurs .== 0;

    llon, llat = ndgrid(longrid, latgrid)
    myproj = ccrs.PlateCarree()
    fig = PyPlot.figure(figsize=(12,8))
    ax = PyPlot.subplot(111, projection=myproj)
    ax.plot(lonobs[data_presence], latobs[data_presence], "wo", markersize=1., zorder=3)
    ax.plot(lonobs[data_absence], latobs[data_absence], "ko", markersize=1., zorder=3)
    pcm = ax.pcolor(llon, llat, dens, cmap=PyPlot.cm.hot_r, zorder=2)
    colorbar(pcm, orientation="horizontal"  )
    decorate_map(ax)
    title(titletext)
    if length(figname) > 0
        PyPlot.savefig(figname, dpi=300, bbox_inches="tight")
        PyPlot.close()
    else
        PyPlot.show()
    end
end
