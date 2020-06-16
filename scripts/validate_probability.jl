using Interpolations

function interpcv(xyi,value_analysis,xy)
    ncv = length(xy[1])
    value_analysis_cv = zeros(ncv)

    xyi = map(x -> collect(x),xyi);

    tmp_itp = interpolate(xyi, value_analysis, Gridded(Linear()))
    itp = extrapolate(tmp_itp,Line())

    #save("tmp.jld2","xyi",xyi,"value_analysis",value_analysis,"xy",xy)
    xyi = zeros(length(xy))

    for i = 1:ncv
    	for j = 1:length(xy)
	        xyi[j] = xy[j][i]
        end

        value_analysis_cv[i] = itp(xyi...)
    end

    return value_analysis_cv
end

"""
    validate_probability(xyi,analysis_probability,xy_cv,occurs_cv)

The function returns the  log-likelihood for a probabilistic model for binary
classification. `xyi` is a tuple with the analysis grid (e.g. `gridlon`,
`gridlat`), `analysis_probability` is the gridded probability of occurance
`xy_cv` is a tuple with the cross-validation points `(obslon,obslat)` and
`occurs_cv` a vector of zeros (absence) and ones (presence).
"""
function validate_probability(xyi,analysis_probability,xy_cv,occurs_cv )
    if isempty(occurs_cv)
        return NaN
    end

    #@show extrema(analysis_probability[.!isnan.(analysis_probability)])
    p = interpcv(xyi,analysis_probability,xy_cv)

    sels = .!isnan.(p)

    y = occurs_cv[sels]
    p = p[sels]

    J = -sum(y .* log.(p)  + (1 .- y) .* log.(1 .- p)) / length(y)
    return J
end
