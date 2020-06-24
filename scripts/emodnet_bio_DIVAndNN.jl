using DIVAnd
using DIVAndNN
using LinearAlgebra
using Statistics
using Random

include("emodnet_bio_grid.jl")
include("validate_probability.jl")
include("PhytoInterp.jl")

if !isdir(datadir)
    mkpath(datadir)
    mkpath(joinpath(datadir,"tmp"))
end

datafile = joinpath(datadir, "Biddulphia sinensis-1995-2020.csv")

# download

function maybedownload(url,fname)
    if !isfile(fname)
        cp(download(url),fname)
    else
        @info("$url is already downloaded")
    end
end

maybedownload("https://dox.ulg.ac.be/index.php/s/VgLglubaTLetHzc/download",
              joinpath(datadir,"gebco_30sec_4.nc"))

maybedownload("https://ec.oceanbrowser.net/data/emodnet-projects/Phase-3/Combined/Water_body_phosphate_combined_V1.nc",
              joinpath(datadir,"tmp","Water_body_phosphate_combined_V1.nc"))

maybedownload("https://ec.oceanbrowser.net/data/emodnet-projects/Phase-3/Combined/Water_body_nitrogen_combined_V1.nc",
              joinpath(datadir,"tmp","Water_body_nitrogen_combined_V1.nc"))

maybedownload("https://ec.oceanbrowser.net/data/emodnet-projects/Phase-3/Combined/Water_body_silicate_combined_V1.nc",
              joinpath(datadir,"tmp","Water_body_silicate_combined_V1.nc"))

maybedownload("https://dox.ulg.ac.be/index.php/s/VgLglubaTLetHzc/download", datafile)

bathname = joinpath(datadir,"gebco_30sec_4.nc");
maskname = joinpath(datadir,"mask.nc")
bathisglobal = true

DIVAndNN.prep_mask(bathname,bathisglobal,gridlon,gridlat,years,maskname)
DIVAndNN.prep_tempsalt(gridlon,gridlat,data_TS,datadir)
DIVAndNN.prep_bath(bathname,bathisglobal,gridlon,gridlat,datadir)
 

BLAS.set_num_threads(1)


ds = Dataset(maskname,"r")
mask = nomissing(ds["mask"][:,:]) .== 1
close(ds)

bathisglobal = true

mask_unused,pmn,xyi = DIVAnd.domain(bathname,bathisglobal,gridlon,gridlat);

# load covariables
covars_fname = [
    ("bathymetry.nc","batymetry",identity),
    ("nitrogen.nc",      "nitrogen",identity),
    ("phosphate.nc",     "phosphate",identity),
    ("silicate.nc",      "silicate",identity),
]

ncovars = length(covars_fname) + 1

sz = (length(gridlon),length(gridlat),ncovars)

field = zeros(sz)

for i = 1:length(covars_fname)
    (fname,varname,trans) = covars_fname[i]

    Dataset(joinpath(datadir,fname)) do ds
        tmp = nomissing(ds[varname][:],NaN)
        tmp = trans.(tmp)
        if ndimensions == 3
            if ndims(tmp) == 2
                field[:,:,:,i] = repeat(tmp,inner = (1,1,length(years)))
            else
                field[:,:,:,i] = tmp
            end
        else
            field[:,:,i] = tmp
        end
    end
end


X,Y = DIVAnd.ndgrid(gridlon,gridlat)

field[:,:,end] .= 1


@show size(field)
function std_or_1(tmp)
    s = std(tmp)
    if s == 0
        return one(s)
    else
        return s
    end
end

# normalize

for n = 1:size(field,3)
    tmp = field[:,:,n][mask];
    field[:,:,n] = (field[:,:,n] .- mean(tmp)) ./ std_or_1(tmp)
end


data_analysis = DIVAndNN.Format2020(datadir,"")

scientificname_accepted = listnames(data_analysis);

lent = 0.6 # years
lent = 0. # years
niter = 500
trainfrac = 0.01
epsilon2ap = 10
epsilon2_background = 10
NLayers = [size(field)[end],4,1]
learning_rate = 0.001
L2reg = 0.0001
dropoutprob = 0.6
len = 75e3


outdir = joinpath(datadir,"Results","emodnet-bio-2020-ncovars$(length(covars_fname))-epsilon2ap$(epsilon2ap)-len$(len)-niter$(niter)-nlayers$(length(NLayers))")
mkpath(outdir)

sname = String(scientificname_accepted[1])

@info sname

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
