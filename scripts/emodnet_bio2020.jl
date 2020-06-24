#using PyPlot
using CSV
using DataFrames
using Dates
using DelimitedFiles
import DIVAnd
using FileIO
using Glob
using Interpolations
using JSON
using Missings
using NCDatasets
using Printf
using Proj4
using Random
using Statistics
using Base.Threads
using LinearAlgebra
using DIVAndNN


BLAS.set_num_threads(1)


include("validate_probability.jl")
include("PhytoInterp.jl")
include("emodnet_bio_grid.jl")

Random.seed!(1234)

# land-sea mask and domain parameters

maskname = joinpath(datadir,"mask.nc");

ds = Dataset(maskname,"r")
mask = nomissing(ds["mask"][:,:]) .== 1
close(ds)

if ndimensions == 3
    mask = repeat(mask,inner=(1,1,length(years)))
end

bathname = joinpath(datadir,"gebco_30sec_4.nc");
bathisglobal = true;


if ndimensions == 3
    mask2,pmn,xyi = DIVAnd.domain(bathname,bathisglobal,gridlon,gridlat,years);
else
    mask2,pmn,xyi = DIVAnd.domain(bathname,bathisglobal,gridlon,gridlat);
end


covars_coord = false
covars_const = true

# load covariables
covars_fname = [
    ("bathymetry.nc","batymetry",identity),
    #            ("dist2coast_subset.nc","distance",identity),
    #("Chlorophyll/chloro_reinterp.nc","chla",identity),
    #("oxygen_reinterp2.nc","oxygen",identity),
    #                ("salinity.nc","salinity",log),
    #               ("temperature.nc","temperature",identity),
    ("nitrogen.nc",      "nitrogen",identity),
    ("phosphate.nc",     "phosphate",identity),
    ("silicate.nc",      "silicate",identity),
]
#covars_fname = []

ncovars = length(covars_fname)

if covars_const
    ncovars = ncovars + 1
end

if covars_coord
    ncovars = ncovars + ndimensions
end

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

if covars_const
    i = length(covars_fname)+1
    @info "add const covariable"

    if ndimensions == 3
        field[:,:,:,i] .= 1
    else
        field[:,:,i] .= 1
    end
end

if covars_coord
    if ndimensions == 3
        @info "add lon/lat/time as covariable"
        field[:,:,:,end-2] = repeat(X,inner = (1,1,length(years)))
        field[:,:,:,end-1] = repeat(Y,inner = (1,1,length(years)))
        field[:,:,:,end]   = repeat(reshape(years,(1,1,length(years))),inner = (length(gridlon),length(gridlat),1))
    else
        @info "add lon/lat as covariable"
        field[:,:,end-1] = X
        field[:,:,end] = Y
    end
end

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

for n = 1:size(field,ndimensions+1)
    if ndimensions == 3
        field[:,:,:,n] = (field[:,:,:,n] .- mean(tmp)) ./ std_or_1(tmp)
    else
        tmp = field[:,:,n][mask];
        field[:,:,n] = (field[:,:,n] .- mean(tmp)) ./ std_or_1(tmp)
    end
end


data_analysis = DIVAndNN.Format2020(expanduser("~/tmp/Emodnet-Bio2020/CSV-split"),"analysis")
data_validation = DIVAndNN.Format2020(expanduser("~/tmp/Emodnet-Bio2020/CSV-split"),"validation")


scientificname_accepted = listnames(data_analysis);

lent = 0.6 # years
lent = 0. # years
niter = 100000
#niter = 100000
#niter = 300000
#niter = 10
niter = 2000*100
#niter = 400000 * 16
#testing
#niter = 10
#niter = 2000
#niter = 2000*20
niter = 500
trainfrac = 0.01
#for l = 1:Ntries
l=1


#epsilon2ap = 1.5
#epsilon2ap = 2
epsilon2ap = 5
epsilon2ap = 10
#epsilon2ap = 50
#@show std(value)
#epsilon2ap = epsilon2
#epsilon2ap = 0.5
#epsilon2ap = 1.
#epsilon2ap = 0.01
#epsilon2ap = 0.1
#epsilon2ap = 0.001
#epsilon2ap = 0.005


epsilon2_background = 10

#NLayers = [size(field)[end],3,1]
NLayers = [size(field)[end],4,1]

#NLayers = [size(field)[end],4,2,1]
#NLayers = [size(field)[end],5,1]
#NLayers = [size(field)[end],2,1]
#NLayers = [size(field)[end],1]
#NLayers = []

learning_rate = 0.00001
learning_rate = 0.001
L2reg = 0.0001
dropoutprob = 0.01
dropoutprob = 0.1
dropoutprob = 0.6
#dropoutprob = 0.99

len = 50e3
len = 50e3
len = 80e3
len = 75e3
#len = 100e3
#len = 120e3
#len = 150e3
#len = 30e3
#len = 20e3

#for len = [50e3, 75e3, 100e3, 125e3]
#    for epsilon2ap = [1, 5, 10, 50, 100]

        outdir = joinpath(datadir,"Results","emodnet-bio-2020-ncovars$(length(covars_fname))-epsilon2ap$(epsilon2ap)-len$(len)-niter$(niter)-nlayers$(length(NLayers))")
        mkpath(outdir)

        nameindex = parse(Int,get(ENV,"INDEX","1"))

        #Threads.@threads for nameindex in 1:length(scientificname_accepted)
        for nameindex in 1:length(scientificname_accepted)

            sname = String(scientificname_accepted[nameindex])
            #sname = "Lithodesmium undulatum"
#sname = "Asterionella kariana"
#sname = "Chaetoceros decipiens"
            global loss_iter
            global val_iter
            @info sname

            paramname = joinpath(outdir,"DIVAndNN_$(sname)_interp.json")

            #if isfile(paramname)
            #    continue
            #end

            lon_a,lat_a,obstime_a,value_a,ids_a = loadbyname(data_analysis,years,sname)
            lon_cv,lat_cv,obstime_cv,value_cv,ids_cv = loadbyname(data_validation,years,sname)

            #time = Float64.(Dates.year.(obstime_a))

            @show value_a[1:min(end,10)]
            @show length(value_a)

            Random.seed!(1234)

            value_analysis = zeros(size(mask))

            xobs_a = if ndimensions == 3
                (lon_a,lat_a,time_a)
            else
                (lon_a,lat_a)
            end

            lenxy = if ndimensions == 3
                (len,len,lent)
            else
                (len,len)
            end

            loss_iter = []
            val_iter = []

            #value_a[(lon_a .< 4) .& (lat_a .< 52)] .= 1.
            #value_a .= 1
            function plotres(i,lossi,value_analysis,y,gradloss,out,iobssel,obspos)
                vp = validate_probability((gridlon,gridlat),value_analysis,(lon_cv,lat_cv),value_cv)
                push!(loss_iter,lossi)
                push!(val_iter,vp)
	            @printf("| %10d | %30.5f | %30.5f |\n",i,lossi,vp)
            end

            value_analysis[:],fw0 = DIVAndNN.analysisprob(
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
	            plotres = plotres,
	            plotevery = 100,
                rmaverage = true,
                trainfrac = trainfrac,
                epsilon2_background = epsilon2_background,
            )

            vp = validate_probability((gridlon,gridlat),value_analysis,(lon_cv,lat_cv),value_cv)
            @show vp

            outname = joinpath(outdir,"DIVAndNN_$(sname)_interp.nc")

            create_nc_results(outname, gridlon, gridlat, value_analysis, sname;
                              varname = "probability", long_name="occurance probability");

            open(paramname,"w") do f
                write(f,JSON.json(
                    Dict(
                        "validation" => vp,
                        "L2reg" =>            L2reg,
                        "dropoutprob" =>      dropoutprob,
                        "epsilon2ap" =>       epsilon2ap,
                        "epsilon2_background" => epsilon2_background,
                        "len" =>              len,
                        "niter" =>            niter,
                        "learning_rate" =>    learning_rate,
                        "NLayers" =>    NLayers,
                        "name" =>    sname,
                        "loss_iter" => loss_iter,
                        "val_iter" => val_iter,
                        "covars" => first.(covars_fname),
                    )
                ))
            end

        end

        score = DIVAndNN.summary(outdir)

        paramname2 = joinpath(outdir,"DIVAndNN.json")

        open(paramname2,"w") do f
            write(f,JSON.json(
                Dict(
                    "validation" => score,
                    "L2reg" =>            L2reg,
                    "dropoutprob" =>      dropoutprob,
                    "epsilon2ap" =>       epsilon2ap,
                    "epsilon2_background" => epsilon2_background,
                    "len" =>              len,
                    "niter" =>            niter,
                    "learning_rate" =>    learning_rate,
                    "NLayers" =>    NLayers,
                    "covars" => first.(covars_fname),
                )
            ))
        end

#    end
#end
