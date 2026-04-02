#using Statistics

function omitnan(x::AbstractVector)
    return filter(!isnan, x)
end

function funcf(f::Function, data::AbstractVector, nan_opt::Int = 1)
    if isempty(data)
        return NaN
    else
        if nan_opt == 1
            data_tmp = omitnan(data)
            return isempty(data_tmp) ? NaN : f(data_tmp)
        else
            return f(data)
        end
    end
end

meanf(data::AbstractVector) = funcf(mean, data)
stdf(data::AbstractVector) = funcf(std, data)
varf(data::AbstractVector) = funcf(var, data)
sumf(data::AbstractVector) = funcf(sum, data)
function covf(
    x::AbstractVector, 
    y::AbstractVector, 
    nan_opt::Int = 1
)
    if length(x) != length(y)
        @warn "len_err: xlen = $(length(x)), ylen = $(length(y))"
        return NaN
    end
    if isempty(x) || isempty(y)
        return NaN
    else
        if nan_opt == 1
            nan_x_flag = isnan.(x)
            nan_y_flag = isnan.(y)
            nan_flag = nan_x_flag .| nan_y_flag
            if all(nan_flag)
                return NaN
            else
                x = x[.!nan_flag]
                y = y[.!nan_flag]
                return cov(x, y)
            end
        else
            return cov(x, y)
        end
    end
end
