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
    plot_bathymetry(bx, by, b, figname)

Plot the bathymetry and save it as `figname`.

## Examples
```julia-repl
julia> plot_bathymetry(bx, by, b, "northsea_bathymetry.png")
```
"""
function plot_bathymetry(bx, by, b, figname="";
    domain=[-16., 9., 45., 66.], dlon=2., dlat=2.)

    fig = PyPlot.figure()
    ax = PyPlot.subplot(111, projection=myproj)
    pcm = pcolor(bx, by, b', vmin=0.);
    colorbar(pcm, orientation="vertical")
    decorate_map(ax, true; domain=domain, dlon=dlon, dlat=dlat)
    title("Depth")
    if length(figname) > 0
        PyPlot.savefig(figname, dpi=300, bbox_inches="tight")
    end
    PyPlot.show()
    PyPlot.close()
end

"""
    plot_mask(bx, by, mask, figname)

Plot the bathymetry and save it as `figname`.

## Examples
```julia-repl
julia> plot_mask(bx, by, mask, "northsea_mask.png")
```
"""
function plot_mask(bx, by, mask, figname="";
    domain=[-16., 9., 45., 66.], dlon=2., dlat=2.)

    fig = PyPlot.figure()
    ax = PyPlot.subplot(111, projection=myproj)
    pcolor(bx, by, mask', cmap=PyPlot.cm.binary_r);
    title("Land-sea mask")
    decorate_map(ax, false; domain=domain, dlon=dlon, dlat=dlat)
    if length(figname) > 0
        PyPlot.savefig(figname, dpi=300, bbox_inches="tight")
    end
    PyPlot.show()
    PyPlot.close()
end



"""
    make_plot_presence_absence(lon, lat, figname)

Plot the locations of absence and presence in 2 subplots and save it as `figname`.

## Examples
```julia-repl
julia> make_plot_presence_absence(lon, lat, "data_locations.png")
```
"""
function make_plot_presence_absence(lon::Vector, lat::Vector, occurs, figname::String="",
    domain=[-16., 9., 45., 66.], dlon=2., dlat=2.)

    data_presence = occurs .== 1;
    data_absence = .!(data_presence);

    fig = PyPlot.figure(figsize=(8, 8))
    ax = PyPlot.subplot(211, projection=myproj)
    ax.plot(lon[data_presence], lat[data_presence], "ko", markersize=.2)
    title("Presence")
    decorate_map(ax, true; domain=domain, dlon=dlon, dlat=dlat)

    ax = PyPlot.subplot(212, projection=myproj)
    ax.plot(lon[data_absence], lat[data_absence], "ko", markersize=.2)
    title("Absence")
    decorate_map(ax, true; domain=domain, dlon=dlon, dlat=dlat)
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
function decorate_map(ax, plotcoast=true;
         domain=[-16., 9., 45., 66.], dlon::Float64=2., dlat::Float64=2.)
    PyPlot.grid(linewidth=0.2)
    if plotcoast
        ax.add_feature(coast, color=".6",
            edgecolor="k", zorder=5)
    end
    ax.set_xlim(domain[1], domain[2])
    ax.set_ylim(domain[3], domain[4])
    ax.set_xticks(domain[1]:dlon:domain[2])
    ax.set_yticks(domain[3]:dlat:domain[4])
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
end

"""
    plot_heatmap(longrid, latgrid, dens, lonobs, latobs, occurs, titletext, figname;
                 domain, dlon, dlat)

Plot the heatmap `dens` (computed by `DIVAnd`) on the grid defined by `longrid`, `latgrid`
and add the locations of the observations defined by `lonobs`, `latobs` and `occurs`.

The vector `occurs` (0's and 1's) is used to separate 'presence' and 'absence' data points.

## Examples
```julia-repl
julia> plot_heatmap(longrid, latgrid, dens, lonobs, latobs, occurs,
    "Actinocyclus", "Actinocyclus_field.png")
```
"""
function plot_heatmap(longrid::StepRangeLen, latgrid::StepRangeLen,
    dens::Array, lonobs::Vector, latobs::Vector, occurs::Vector,
    titletext::String, figname::String="", vmin=0., vmax=maximum(dens);
    domain=[-16., 9., 45., 66.], dlon=2., dlat=2.)

    data_presence = occurs .== 1;
    data_absence = occurs .== 0;

    llon, llat = ndgrid(longrid, latgrid)
    myproj = ccrs.PlateCarree()
    fig = PyPlot.figure(figsize=(12,8))
    ax = PyPlot.subplot(111, projection=myproj)
    decorate_map(ax, true; domain=domain, dlon=dlon, dlat=dlat)
    ax.plot(lonobs[data_presence], latobs[data_presence], "go", markersize=1, zorder=4,
            label="Presence")
    ax.plot(lonobs[data_absence], latobs[data_absence], "ko", markersize=2.5, zorder=3,
            label="Absence")
    pcm = ax.pcolor(llon, llat, dens, cmap=PyPlot.cm.RdYlBu_r, zorder=2, vmin=vmin, vmax=vmax)
    PyPlot.legend(loc=2)
    colorbar(pcm, orientation="vertical")
    title(titletext)
    if length(figname) > 0
        PyPlot.savefig(figname, dpi=300, bbox_inches="tight")
        PyPlot.close()
    else
        PyPlot.show()
    end
end

"""
    plot_error(longrid, latgrid, error, titletext, figname;
                 domain, dlon, dlat)

Plot the error field `error` (computed by `DIVAnd`) on the grid defined by `longrid`, `latgrid`.

## Examples
```julia-repl
julia> plot_heatmap(longrid, latgrid, error, "Actinocyclus", "Actinocyclus_error.png")
```
"""
function plot_error(longrid::StepRangeLen, latgrid::StepRangeLen,
    error::Array, titletext::String="", figname::String="";
    domain=[-16., 9., 45., 66.], dlon::Float64=2., dlat::Float64=2.)

    llon, llat = ndgrid(longrid, latgrid)
    myproj = ccrs.PlateCarree()
    fig = PyPlot.figure(figsize=(12,8))
    ax = PyPlot.subplot(111, projection=myproj)
    pcm = ax.pcolor(llon, llat, error, zorder=2, cmap=PyPlot.cm.RdYlGn_r)
    colorbar(pcm, orientation="vertical")
    decorate_map(ax, true; domain=domain, dlon=dlon, dlat=dlat)
    title(titletext)
    if length(figname) > 0
        PyPlot.savefig(figname, dpi=300, bbox_inches="tight")
        PyPlot.close()
    else
        PyPlot.show()
    end
end


"""
plot_bincounts(longrid, latgrid, mask, count, speciesname, figname)

Plot the number of observations per grid cell.
The grid is defined by longrid, latgrid.

## Examples
```julia-repl
julia> plot_bincounts(longrid, latgrid, mask, count, "Actinocyclus", "Actinocyclus_bins.png")
```
"""
function plot_bincounts(longrid, latgrid, mask, count::Array,
        speciesname::String="", figname::String="";
        domain=[-16., 9., 45., 66.], dlon=2., dlat=2., vmin=0., vmax=2.)

    xx = [i for i in longrid, j in latgrid]
    yy = [j for i in longrid, j in latgrid];

    PyPlot.figure(figsize=(10, 8))
    ax = PyPlot.subplot(111, projection=myproj)
    ax.pcolor(xx, yy, mask, cmap=PyPlot.cm.binary_r)
    pcm = ax.pcolor(xx, yy, count, cmap=PyPlot.cm.hot_r, vmin=vmin, vmax=vmax)
    cb = colorbar(pcm, extend="max", shrink=.8)

    if length(speciesname) > 0
        PyPlot.title("Logarithm of the number observations per bin for $(speciesname)")
    end

    decorate_map(ax, false; domain=domain, dlon=dlon, dlat=dlat)
    if length(figname) > 0
        savefig(figname, dpi=300, bbox_inches="tight")
        PyPlot.close()
    else
        PyPlot.show()
    end

end
