using PyCall
using PyPlot
using Proj4
using DIVAnd
using NCDatasets

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

function read_data_phyto(datafile::String)
    data = readdlm(datafile, ',');
    columntitles = data[1,:]

    df = DateFormat("y-m-d");
    dates = Date.(data[2:end,6], df)
    lons = Float64.(data[2:end,7])
    lats = Float64.(data[2:end,8])
    occurs = data[2:end,end-3];

    return dates::Vector{Date}, lons::Vector{Float64}, lats::Vector{Float64}, occurs::Vector
end

"""
    get_species_name(datafile)

Return the species name based on the data file name
"""
function get_species_name(datafile::String)::String
    return split(datafile, "-")[1]
end

function transform_coords(lon::Array, lat::Array)

    # Setup projections
    wgs84 = Projection("+proj=longlat +datum=WGS84 +no_defs")
    espgs32361 = Projection("+proj=utm +zone=31 +north +datum=WGS84 +units=m +no_defs")

    # Allocate
    npoints = length(lons)
    lonp = Array{Float64, 1}(undef, npoints)
    latp = Array{Float64, 1}(undef, npoints)

    # Loop
    for i = 1:npoints
        lonp[i], latp[i], e = transform(espgs32361, wgs84, [lon[i], lat[i], 0.])
    end

    return lonp, latp
end

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

function create_nc_results(filename::String, lons, lats, field,
    speciesname::String=""; valex=-999.9)
    Dataset(filename, "c") do ds

        # Dimensions
        ds.dim["lon"] = length(lons)
        ds.dim["lat"] = length(lats)
        #ds.dim["time"] = Inf # unlimited dimension

        # Declare variables
        ncfield = defVar(ds,"heatmap", Float64, ("lon", "lat"))
        ncfield.attrib["missing_value"] = Float64(valex)
        ncfield.attrib["_FillValue"] = Float64(valex)
        ncfield.attrib["long_name"] = "Interpolated abundance"
        ncfield.attrib["units"] = "ind/m3"

        """
        nctime = defVar(ds,"time", Float32, ("time",))
        nctime.attrib["missing_value"] = Float32(valex)
        nctime.attrib["units"] = "seconds since 1981-01-01 00:00:00"
        nctime.attrib["time"] = "time"
        """

        nclon = defVar(ds,"lon", Float32, ("lon",))
        nclon.attrib["missing_value"] = Float32(valex)
        nclon.attrib["_FillValue"] = Float32(valex)
        nclon.attrib["units"] = "degrees East"
        nclon.attrib["lon"] = "longitude"

        nclat = defVar(ds,"lat", Float32, ("lat",))
        nclat.attrib["missing_value"] = Float32(valex)
        nclat.attrib["_FillValue"] = Float32(valex)
        nclat.attrib["units"] = "degrees North"
        nclat.attrib["lat"] = "latitude"

        # Global attributes
        ds.attrib["institution"] = "GHER - University of Liege"
        ds.attrib["title"] = "Heatmap based on abundance of $(speciesname)"
        ds.attrib["comment"] = "Original data prepared by Deltares"
        ds.attrib["data authors"] = "Luuk van der Heijden (Luuk.vanderHeijden@deltares.nl), Willem Stolte (Willem.Stolte@deltares.nl)"
        ds.attrib["processing authors"] = "C. Troupin (ctroupin@uliege), A. Barth (a.barth@uliege.be)"
        ds.attrib["created"] = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")

        # Define variables
        ncfield[:] = field
        nclon[:] = lons
        nclat[:] = lats;

    end
end;

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


function plot_heatmap(longrid::StepRangeLen, latgrid::StepRangeLen,
    dens::Array, lonobs::Vector, latobs::Vector,
    titletext::String, figname::String="")

    llon, llat = ndgrid(longrid, latgrid)
    myproj = ccrs.PlateCarree()
    fig = PyPlot.figure(figsize=(12,8))
    ax = PyPlot.subplot(111, projection=myproj)
    ax.plot(lonobs, latobs, "ko", markersize=0.2)
    pcm = ax.pcolor(llon, llat, dens, cmap=PyPlot.cm.hot_r)
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
