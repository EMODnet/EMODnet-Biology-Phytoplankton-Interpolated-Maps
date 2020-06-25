using DIVAnd
using DIVAndNN
using LinearAlgebra
using Statistics
using Random
using Dates

include("emodnet_bio_grid.jl")
include("validate_probability.jl")
include("PhytoInterp.jl")

# create directories
if !isdir(datadir)
    mkpath(datadir)
    mkpath(joinpath(datadir,"tmp"))
end


# download file from URL is necessary

function maybedownload(url,fname)
    if !isfile(fname)
        cp(download(url),fname)
    else
        @info("$url is already downloaded")
    end
end


# GEBCO Bathymetry
bathname = joinpath(datadir,"gebco_30sec_4.nc");
maybedownload("https://dox.ulg.ac.be/index.php/s/VgLglubaTLetHzc/download",
              joinpath(datadir,"gebco_30sec_4.nc"))


# Sample data file
datafile = joinpath(datadir, "Biddulphia sinensis-1995-2020.csv")
maybedownload("https://dox.ulg.ac.be/index.php/s/VgLglubaTLetHzc/download", datafile)

maskname = joinpath(datadir,"mask.nc")
bathisglobal = true


# Environmental covariables

#=
maybedownload("https://ec.oceanbrowser.net/data/emodnet-projects/Phase-3/Combined/Water_body_phosphate_combined_V1.nc",
              joinpath(datadir,"tmp","Water_body_phosphate_combined_V1.nc"))

maybedownload("https://ec.oceanbrowser.net/data/emodnet-projects/Phase-3/Combined/Water_body_nitrogen_combined_V1.nc",
              joinpath(datadir,"tmp","Water_body_nitrogen_combined_V1.nc"))

maybedownload("https://ec.oceanbrowser.net/data/emodnet-projects/Phase-3/Combined/Water_body_silicate_combined_V1.nc",
              joinpath(datadir,"tmp","Water_body_silicate_combined_V1.nc"))

DIVAndNN.prep_mask(bathname,bathisglobal,gridlon,gridlat,years,maskname)
DIVAndNN.prep_tempsalt(gridlon,gridlat,data_TS,datadir)
DIVAndNN.prep_bath(bathname,bathisglobal,gridlon,gridlat,datadir)

=#

maybedownload("https://dox.ulg.ac.be/index.php/s/y9Z0c1wb5YshVDW/download",
              joinpath(datadir,"silicate.nc"))

maybedownload("https://dox.ulg.ac.be/index.php/s/A1NPSWwQYkx6Wy6/download",
              joinpath(datadir,"phosphate.nc"))

maybedownload("https://dox.ulg.ac.be/index.php/s/LDPbPWBvW6wPmCw/download",
              joinpath(datadir,"nitrogen.nc"))


BLAS.set_num_threads(1)


# load mask
ds = Dataset(maskname,"r")
mask = nomissing(ds["mask"][:,:]) .== 1
close(ds)


# compute local resolution
mask_unused,pmn,xyi = DIVAnd.domain(bathname,bathisglobal,gridlon,gridlat);

# load covariables
covars_fname = [
    # filename       ,  var. name  , transformation function
    ("bathymetry.nc" , "batymetry" , identity),
    ("nitrogen.nc"   , "nitrogen"  , identity),
    ("phosphate.nc"  , "phosphate" , identity),
    ("silicate.nc"   , "silicate"  , identity),
]

covars_fname = map(entry -> (joinpath(datadir,entry[1]),entry[2:end]...),covars_fname)

field = DIVAndNN.loadcovar((gridlon,gridlat),covars_fname;
                           covars_const = true)

# normalize covariables
DIVAndNN.normalize!(mask,field)


# Inventory of all data files
data_analysis = DIVAndNN.Format2020(datadir,"")
scientificname_accepted = listnames(data_analysis);

# Parameters
niter = 500
trainfrac = 0.01
epsilon2ap = 10
epsilon2_background = 10
NLayers = [size(field)[end],4,1]
learning_rate = 0.001
L2reg = 0.0001
dropoutprob = 0.6
len = 75e3

# output directory
outdir = joinpath(datadir,"Results","emodnet-bio-2020-ncovars$(length(covars_fname))-epsilon2ap$(epsilon2ap)-len$(len)-niter$(niter)-nlayers$(length(NLayers))")
mkpath(outdir)

sname = String(scientificname_accepted[1])

@info sname

# load variable
lon_a,lat_a,obstime_a,value_a,ids_a = loadbyname(data_analysis,years,sname)

Random.seed!(1234)

xobs_a = (lon_a,lat_a)
lenxy = (len,len)

value_analysis,fw0 = DIVAndNN.analysisprob(
    mask,pmn,xyi,xobs_a,
    value_a,
    lenxy,epsilon2ap,
    field,
    NLayers,
    costfun = DIVAndNN.nll,
    niter = niter,
    dropoutprob = dropoutprob,
    L2reg = L2reg,
    learning_rate = learning_rate,
	plotevery = 100,
    rmaverage = true,
    trainfrac = trainfrac,
    epsilon2_background = epsilon2_background,
)

outname = joinpath(outdir,"DIVAndNN_$(sname)_interp.nc")

create_nc_results(outname, gridlon, gridlat, value_analysis, sname;
                  varname = "probability", long_name="occurance probability");


include("emodnet_bio_plot2.jl")
