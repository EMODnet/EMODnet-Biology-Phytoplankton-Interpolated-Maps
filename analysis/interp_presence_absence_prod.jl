using DIVAnd
using PyPlot
using Proj4
using DelimitedFiles
using PyCall
using Dates
using NCDatasets
include("../scripts/PhytoInterp.jl")
include("../scripts/PhytoInterpPlot.jl")
doplot = true
writenc = false 

if doplot
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
end 

# csvdir = "/data/EMODnet/Biology/phytoplankton/csv/" # North Sea 
csvdir = "/media/ctroupin/My Passport/data/EMODnet/Biology/phytoplankton/GreaterNorthSea/csv_files/" #Greater North Sea
@info(isdir(csvdir));

datadir = "../data/"
datafilelist = readdir(csvdir);
@info("Working on $(length(datafilelist)) files");
figdir = "../product/figures/heatmap/greaterNorthSea/" 
binfigdir = "../product/figures/bins/greaterNorthSea/"
datafigdir = "../product/figures/data/greaterNorthSea/"
outputdir = "../product/netCDF/greaterNorthSea/Phytoplankton/V1/"
isdir(datadir) ? " " : mkpath(datadir)
isdir(datafigdir) ? " " : mkpath(datafigdir)
isdir(figdir) ? " " : mkpath(figdir);
isdir(binfigdir) ? " " : mkpath(binfigdir);
isdir(outputdir) ? " " : mkpath(outputdir);

# North Sea
#longrid = -2.:0.1:10.
#latgrid = 51.:0.1:56.

# Greater North Sea
longrid = -16.:0.1:9.
latgrid = 45.:0.1:66.

bathname = joinpath(datadir, "gebco_30sec_4.nc")
if !isfile(bathname)
    download("https://dox.ulg.ac.be/index.php/s/RSwm4HPHImdZoQP/download", bathname)
else
    @info("Bathymetry file already downloaded")
end

bx, by, b = load_bath(bathname, true, longrid, latgrid)
@show size(b)
if doplot
    plot_bathymetry(bx, by, b, joinpath(figdir, "northsea_bathy");
                    domain=[-16., 9., 45., 66.], dlon=5., dlat=5.)
end

mask, (pm, pn),(xi, yi) = DIVAnd.DIVAnd_rectdom(longrid, latgrid);
xi, yi, mask = DIVAnd.load_mask(bathname, true, longrid, latgrid, 0.0);
xx, yy = ndgrid(xi, yi);
if doplot
    plot_mask(bx, by, mask, joinpath(figdir, "northsea_mask"); 
              domain=[-16., 9., 45., 66.], dlon=5., dlat=3.)
end

include("../scripts/PhytoInterpPlot.jl")

for datafile in datafilelist[21:40]
    speciesname = get_species_name(basename(datafile))
    speciesslug = get_species_slug(basename(datafile))
    @info("Working on $(speciesname)")
    @info(speciesslug)
        
    # Data reading
    dates, lons, lats, occurs = read_data_phyto(joinpath(csvdir, datafile));
    
    # Conversion of coordinates
    # We have to go from `EPSG:32361` to `WGS84`.         
    lon, lat = transform_coords(lons, lats)
    
    # Plot presence/absence
    if doplot
        @info("Plotting")
        make_plot_presence_absence(lon, lat, occurs, 
            joinpath(datafigdir, "$(speciesslug)_data_presence_absence.png"), 
            [-16., 9., 45., 66.], 4., 2.)
    end
    
    # Create bins
    obscount, obscountlog = count_obs(longrid, latgrid, lon, lat);
    if doplot
        figname = joinpath(binfigdir, "$(speciesslug)_bins.png")
        plot_bincounts(longrid, latgrid, mask, obscountlog, speciesname, figname; 
                       domain=[-16., 9., 45., 66.], dlon=4., dlat=2., vmax=3.)
    end
    
    
    data_presence = occurs .== 1;
    data_absence = .!(data_presence);
    npre = sum(data_presence)
    nabs = sum(data_absence);
    
    @info("Presence data: $(npre), absence data: $(nabs)")
    
    # Compute heatmap
    inflation = ones(length(lon));
    L = 0.5

    
    @time dens1, LHM, LCV, LSCV = DIVAnd_heatmap(mask, (pm,pn), (xx, yy), 
        (lon, lat), inflation, L);
    @time dens2, LHM, LCV, LSCV = DIVAnd_heatmap(mask, (pm,pn), (xx, yy), 
        (lon[data_presence], lat[data_presence]), inflation[data_presence], L);
    @time dens3, LHM, LCV, LSCV = DIVAnd_heatmap(mask, (pm,pn), (xx, yy), 
        (lon[data_absence], lat[data_absence]), inflation[data_absence], L);
    
    reldens = npre .* dens2 ./ (npre .* dens2 .+ nabs .* dens3);
    reldens2 = npre .* dens2 ./ ( (npre + nabs) .* dens1);
    
    # Compute error field with CPME
    cpme = DIVAnd_cpme(mask, (pm, pn), (xx, yy), (lon, lat), ones(length(lon)), 0.2, 10.);
    
    if doplot
        
        
        plot_heatmap(longrid, latgrid, dens2, lon, lat, occurs,
            "$(speciesname): heatmap with 'presence' data points", 
            joinpath(figdir, "$(speciesslug)_heatmap_presence.png"), 0., 0.2,
            domain=[-16., 9., 45., 66.], dlon=4., dlat=2.)
        
        
        plot_heatmap(longrid, latgrid, dens3, lon, lat, occurs,
            "$(speciesname): heatmap with 'absence' data points", 
            joinpath(figdir, "$(speciesslug)_heatmap_absence.png"), 0., 0.2,
            domain=[-16., 9., 45., 66.], dlon=4., dlat=2.)
        
        """
        plot_heatmap(longrid, latgrid, dens3, lon, lat, occurs,
            "$(speciesname): heatmap with all data points", 
            joinpath(figdir, "$(speciesslug)_heatmap_all.png"), 0., 1.,
            domain=[-16., 9., 45., 66.], dlon=4., dlat=2.)
        """
        
        plot_heatmap(longrid, latgrid, reldens, lon, lat, occurs,
            "$(speciesname) probability", 
            joinpath(figdir, "$(speciesslug)_heatmap_relative.png"), 0., 1.,
            domain=[-16., 9., 45., 66.], dlon=4., dlat=2.)
        
        plot_heatmap(longrid, latgrid, reldens, lon, lat, occurs,
            "$(speciesname) probability", 
            joinpath(figdir, "$(speciesslug)_heatmap_relative2.png"), 0., 1.,
            domain=[-16., 9., 45., 66.], dlon=4., dlat=2.)
    end
    
    # Write in netCDF
    
    if writenc
        create_nc_results(joinpath(outputdir, "$(speciesslug)_heatmap.nc"), 
                          longrid, latgrid, reldens, speciesname);
    
        # Add the error field
        write_nc_error(joinpath(outputdir, "$(speciesslug)_heatmap.nc"), cpme);
    end
    
    figname = joinpath(figdir, "$(speciesslug)_cpme.png")
    @info("Saving image as $(figname)")
    
    if doplot
        plot_error(longrid, latgrid, cpme, "$(speciesname)", 
            joinpath(figdir, "$(speciesslug)_cpme.png"),
            domain=[-16., 9., 45., 66.], dlon=4., dlat=2.)
    end    
    
end


