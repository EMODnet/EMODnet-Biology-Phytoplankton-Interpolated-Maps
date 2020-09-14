using DataFrames, CSV
using Random
using Glob


function splitfile(fname,validation_fraction,csvsplitdir)
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

    CSV.write(joinpath(csvsplitdir,replace(basename(fname),".csv" => ".analysis.csv")),df[.!for_validation,:],quotestrings=true)
    CSV.write(joinpath(csvsplitdir,replace(basename(fname),".csv" => ".validation.csv")),df[for_validation,:],quotestrings=true)
end

"""
This script should be run in the directory containting all original CSV files

# example
```julia
validation_fraction = 0.2
csvsplitdir = "../CSV-split"
csvdir = "./"
splitdir(csvdir,validation_fraction,csvsplitdir)
```
"""
function splitdir(csvdir,validation_fraction,csvsplitdir)
    mkpath(csvsplitdir)

    for fname in glob("*csv",csvdir)
        splitfile(fname,validation_fraction,csvsplitdir)
    end
end



