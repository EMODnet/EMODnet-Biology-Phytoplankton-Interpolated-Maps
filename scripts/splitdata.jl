using DataFrames, CSV
using Random
using Glob

# This script should be run in the directory containting all original CSV files

function splitfile(fname,validation_fraction,splitdir)
    @info "split file $fname"

    # repeatable random number per file
    Random.seed!((codepoint.(collect(fname))))

    df = DataFrame(CSV.read(fname));
    posid = collect(zip(df.xUTM,df.yUTM))

    uposid = unique(posid);
    posid_val = shuffle(uposid)[1: round(Int,validation_fraction*end )];

    @show length(posid_val)
    @show length(uposid)

    for_validation = map(posid -> posid in posid_val,zip(df.xUTM,df.yUTM));

    @show sum(for_validation)
    @show length(for_validation)

    CSV.write(joinpath(splitdir,replace(fname,".csv" => ".analysis.csv")),df[.!for_validation,:],quotestrings=true)
    CSV.write(joinpath(splitdir,replace(fname,".csv" => ".validation.csv")),df[for_validation,:],quotestrings=true)
end

validation_fraction = 0.2

splitdir = "../CSV-split"
csvdir = "./"

for fname in glob("*csv",csvdir)
    splitfile(fname,validation_fraction,splitdir)
end
