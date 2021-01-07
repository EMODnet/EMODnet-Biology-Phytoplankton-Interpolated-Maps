using PyPlot
using Proj4
using DIVAnd
using NCDatasets
using DelimitedFiles
using GridInterpolations

"""
    read_data_phyto(datafile)

Read the coordinates, dates and occurences stored in the file `datafile`

## Examples
```julia-repl
julia> dates, lons, lats, occurs = read_data_phyto("Archaeperidinium-1995-2020.csv")
```
"""
function read_data_phyto(datafile::String)
    data = readdlm(datafile, ',');
    columntitles = data[1,:]

    # Identify column containing the date and the coordinates
    datecolumn = findall(columntitles .== "date")[1]
    xcolumn = findall(columntitles .== "xUTM")[1]
    ycolumn = findall(columntitles .== "yUTM")[1]
    # occurscolumn = findall(columntitles .== "occurs")[1]

    # New data files have column called "occurrence" instead of "occurs"
    # Also they have a column "occurrenceid" that we're not interested in
    occurscolumn = findall(occursin.("occur", columntitles) .& .!occursin.("id", columntitles))[1]

    df = DateFormat("y-m-d");
    dates = Date.(data[2:end,datecolumn], df)
    lons = Float64.(data[2:end,xcolumn])
    lats = Float64.(data[2:end,ycolumn])
    occurs = data[2:end,occurscolumn];

    return dates::Vector{Date}, lons::Vector{Float64}, lats::Vector{Float64}, occurs::Vector
end

"""
    get_species_name(datafile)

Return the species name based on the data file name

## Examples
```julia-repl
julia> speciesname = get_species_name("Archaeperidinium-1995-2020.csv")
"Archaeperidinium"
```
"""
function get_species_name(datafile::String)::String
    m = match(r"(.+)-1995-2020.csv", datafile)
    return m[1]
end

"""
    get_species_slug(filename)

Return the a species 'slug' based on the data file name

## Examples
```julia-repl
julia> slug = get_species_slug("Plagiogramma brockmanni var. brockmanni-1995-2020.csv")
"Plagiogramma_brockmanni_var_brockmanni"
```
"""
function get_species_slug(filename::String)::String
    speciesname = get_species_name(filename)
    speciesslug = replace(replace(speciesname, " " => "_"), "." => "")
    return speciesslug
end

"""
    transform_coords(lon, lat)

Transform the initial coordinates (`lon`, `lat`), in 'UTM zone 31N', into WGS84.

## Examples
```julia-repl
julia> lonp, latp = transform_coords([477850.085718031, 560048.806065514], [5756212.06783869, 5746827.70084791])
([2.677667000000006, 3.872167000000008], [51.956167000000036, 51.86900000000005])
```
"""
function transform_coords(lon::Vector, lat::Vector)

    # Setup projections
    wgs84 = Projection("+proj=longlat +datum=WGS84 +no_defs")
    espgs32361 = Projection("+proj=utm +zone=31 +north +datum=WGS84 +units=m +no_defs")

    # Allocate
    npoints = length(lon)
    lonp = Array{Float64, 1}(undef, npoints)
    latp = Array{Float64, 1}(undef, npoints)

    # Loop
    for i = 1:npoints
        lonp[i], latp[i], e = transform(espgs32361, wgs84, [lon[i], lat[i], 0.])
    end

    return lonp, latp
end

"""
    reinterp_field(longrid, latgrid, field, lonobs, latobs)

Re-interpolate the field defined by `longrid`, `latgrid` and `field` onto the
positions defined by `lonobs`, `latobs`

"""
function reinterp_field(longrid, latgrid, field::Array, lonobs::Vector, latobs::Vector)
    grid = RectangleGrid(longrid, latgrid);
    gridData = field[:];

    # Allocate
    fieldinterp = Vector{Float64}(undef, length(lonobs))
    i = 1
    for (loninterp, latinterp) in zip(lonobs, latobs)
        fieldinterp[i] = GridInterpolations.interpolate(grid, gridData, [loninterp, latinterp])
        i += 1
    end
    return fieldinterp::Vector
end


"""
    create_nc_results(filename, lons, lats, field, speciesname)

Write a netCDF file `filename` with the coordinates `lons`, `lats` and the
heatmap `field`. `speciesname` is used for the 'title' attribute of the netCDF.

## Examples
```julia-repl
julia> create_nc_results("Bacteriastrum_interp.nc", lons, lats, field,
    "Bacteriastrum")
```
"""
function create_nc_results(filename::String, lons, lats, field,
                           speciesname::String="";
                           valex=-999.9,
                           varname = "heatmap",
                           long_name = "Heatmap",
                           time = nothing,
                           )
    Dataset(filename, "c") do ds

        # Dimensions
        ds.dim["lon"] = length(lons)
        ds.dim["lat"] = length(lats)

        # Declare variables

        if time == nothing
            ncfield = defVar(ds, varname, Float64, ("lon", "lat"))
        else
            ds.dim["time"] = length(time)
            ncfield = defVar(ds, varname, Float64, ("lon", "lat", "time"))
            nctime = defVar(ds,"time", Float64, ("time",))
            nctime.attrib["units"] = "years"
            nctime.attrib["time"] = "time"
        end

        ncfield.attrib["missing_value"] = Float64(valex)
        ncfield.attrib["_FillValue"] = Float64(valex)
        ncfield.attrib["long_name"] = long_name

        nclon = defVar(ds,"lon", Float32, ("lon",))
        nclon.attrib["units"] = "degrees East"
        nclon.attrib["lon"] = "longitude"

        nclat = defVar(ds,"lat", Float32, ("lat",))
        nclat.attrib["units"] = "degrees North"
        nclat.attrib["lat"] = "latitude"

        # Global attributes
        ds.attrib["project"] = "EMODnet Biology"
        ds.attrib["institution"] = "GHER - University of Liege"
        ds.attrib["title"] = "$(long_name) based on abundance of $(speciesname)"
        ds.attrib["comment"] = "Original data prepared by Deltares"
        ds.attrib["data authors"] = "Luuk van der Heijden (Luuk.vanderHeijden@deltares.nl), Willem Stolte (Willem.Stolte@deltares.nl)"
        ds.attrib["processing authors"] = "C. Troupin (ctroupin@uliege), A. Barth (a.barth@uliege.be)"
        ds.attrib["created"] = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")

        # Define variables
        ncfield[:] = field

        nclon[:] = lons
        nclat[:] = lats;
        if time != nothing
            nctime[:] = time
        end
    end
end;

"""
    write_nc_error(filename, error)

Write the error field in the file `filename` created by `create_nc_results`

## Examples
```julia-repl
julia> write_nc_error(filename, error)
```
"""
function write_nc_error(filename::String, error::Array  ;
    valex=-999.9, varname = "heatmap",
    long_name = "Heatmap")

    Dataset(filename, "a") do ds
        ncerror = defVar(ds, varname * "_error", Float64, ("lon", "lat"))
        ncerror.attrib["missing_value"] = Float64(valex)
        ncerror.attrib["_FillValue"] = Float64(valex)
        ncerror.attrib["long_name"] = long_name * "_error"
        ncerror[:] = error
    end
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
    titletext::String, figname::String="", vmin=0., vmax=1.)

    data_presence = occurs .== 1;
    data_absence = occurs .== 0;

    llon, llat = ndgrid(longrid, latgrid)
    fig = PyPlot.figure(figsize=(12,8))
    ax = PyPlot.subplot(111)
    ax.plot(lonobs[data_presence], latobs[data_presence], "wo", markersize=1., zorder=3)
    ax.plot(lonobs[data_absence], latobs[data_absence], "ko", markersize=1., zorder=3)
    pcm = ax.pcolor(llon, llat, dens, cmap=PyPlot.cm.hot_r, zorder=2, vmin=vmin, vmax=vmax)
    colorbar(pcm, orientation="horizontal"  )
    title(titletext)
    if length(figname) > 0
        PyPlot.savefig(figname, dpi=300, bbox_inches="tight")
        PyPlot.close()
    else
        PyPlot.show()
    end
end

"""
    count_obs(longrid, latgrid, lon, lat, title)

Count the number of observations `obscount` and its logarithm `obscountlog`
made in each cell of the grid defined by longrid, latgrid.

## Examples
```julia-repl
julia> obscount, obscountlog = count_obs(longrid, latgrid, lon, lat)
```
"""
function count_obs(longrid::StepRangeLen, latgrid::StepRangeLen,
                   lon::Array, lat::Array)

    Δlon = longrid[2] - longrid[1]
    Δlat = latgrid[2] - latgrid[1]

    obscount = Array{Float64, 2}(undef, length(longrid) - 1, length(latgrid) - 1)
    obscountlog = Array{Float64, 2}(undef, length(longrid) -1 , length(latgrid) -1 )

    for (i, long) in enumerate(longrid[1:end-1])
        goodlon = (lon .>= long) .& (lon .< long + Δlon);
        for (j, latg) in enumerate(latgrid[1:end-1])
            goodlat = (lat .>= latg) .& (lat .< latg + Δlat);
            goodpos = findall(goodlon .& goodlat);
            np = length(goodpos);
            if np == 0
                obscount[i, j] = NaN
                obscountlog[i, j] = NaN
            else
                obscount[i, j] = np
                obscountlog[i, j] = log10(np)
            end

        end
    end

    return obscount::Array, obscountlog::Array
end
