#=
caution:: "lag" means "x" is late when lag > 0.
=#
#using ..Config
#using ..LagExporter
function xcov(x_in::AbstractVector, 
    y_in::AbstractVector, 
    lag_min::Int64, 
    lag_max::Int64
)

    @assert !(lag_max < lag_min) "caution: lag_min > lag_max: you must input a correct lag_argument."

    lag_series = lag_min : lag_max
    n_lags = length(lag_series)
    cov_series = fill(NaN, n_lags)

    for (idx, i) in enumerate(lag_series)
        if i <= 0
            x_v = @view x_in[1 : end + i]
            y_v = @view y_in[1 - i : end]
            cov_series[idx] = covf(x_v, y_v)
        else
            x_v = @view x_in[1 + i : end]
            y_v = @view y_in[1 : end - i]
            cov_series[idx] = covf(x_v, y_v)
        end
    end
    
    valid_indices = findall(!isnan, cov_series)
    
    if isempty(valid_indices)
        return NaN, 0
    end
    
    _, tmp_idx = findmax(abs, @view cov_series[valid_indices])
    
    max_idx = valid_indices[tmp_idx]

    return cov_series[max_idx], lag_series[max_idx]
    #=lag = lag_series[abs.(cov_series) .== maximum(abs.(cov_series))]
    if length(lag) == 0
        lag = 0
    else
        lag = lag[1]
    end
    cov_max = cov_series[abs.(cov_series) .== maximum(abs.(cov_series))]
    if length(cov_max) == 0
        cov_max = NaN
    else
        cov_max = cov_max[1]
    end
    lag = Int64(lag)
    return [cov_max, lag]=#

end

function cross_cov(
    x_in::AbstractVector, 
    y_in::AbstractVector, 
    lag_offset::Int, 
    lag_min::Int64, 
    lag_max::Union{Int64, Nothing} = nothing
)
    #the unit of lag_min, lag_max is second * hz.
    if isnothing(lag_max)
        lag_max= -lag_min
    end
    #lag_min = lag_min * SAMPLING_SPAN
    #lag_max = lag_max * SAMPLING_SPAN
    @views begin
        if lag_offset > 0
            x_proc = x_in[1+lag_offset:end]
            y_proc = y_in[1:end-lag_offset]
        else
            x_proc = x_in[1:end+lag_offset]
            y_proc = y_in[1-lag_offset:end]
        end
    end
    res, lag_xy::Int64 = xcov(x_proc, y_proc, lag_min, lag_max)
    return (res, lag_xy::Int64)
end

