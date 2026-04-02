#=
caution:: "lag" means that "x" is late when lag > 0.
=#

#Rannik 2016; Dong et al., 2021
#For calcurating random error for each interval
#Created by Docopy

function lagged_covariances(x_in::Vector{Float64}, y_in::Vector{Float64}, max_shift_sec::Int64, sampling_span::Int64)
    #lag -> max_shift
    @assert max_shift_sec > 100 "caution: too small lag for rannnik2016_rand_err"
    max_shift = max_shift_sec * sampling_span
    min_shift = 100 * sampling_span

    #used shifts
    shifts_neg = -max_shift : -min_shift
    shifts_pos = min_shift : max_shift

    n_neg = length(shifts_neg)
    n_pos = length(shifts_pos)
    covs = Vector{Float64}(undef, n_neg + n_pos)
    len_total = length(x_in)
    #s_x, e_x = max(1, 1 + L), min(len_total, len_total + L)
    #s_y, e_y = max(1 - L, 1), min(len_total - L, len_total)
    #Threads.@threads 
    for i in 1:n_neg
        L = shifts_neg[i]
        x_v = @view x_in[1 : len_total+L]
        y_v = @view y_in[1 - L : len_total]

        n = length(x_v)
        sx = 0.0
        sy = 0.0
        sxy = 0.0
        valid_n = 0
        for j in 1:n
            if !isnan(x_v[j]) && !isnan(y_v[j])
                sx += x_v[j]
                sy += y_v[j]
                sxy += x_v[j] * y_v[j]
                valid_n += 1
            end
        end
        if valid_n > 0
            mx = sx / valid_n
            my = sy / valid_n
            covs[i] = (sxy / valid_n) - (mx * my)
        else
            covs[i] = NaN
        end
    end

    #Threads.@threads
    for i in 1:n_pos
        L = shifts_pos[i]
        x_v = @view x_in[1 + L : len_total]
        y_v = @view y_in[1 : len_total - L]
        n = length(x_v)
        sx = 0.0
        sy = 0.0
        sxy = 0.0
        valid_n = 0
        for j in 1:n
            if !isnan(x_v[j]) && !isnan(y_v[j])
                sx += x_v[j]
                sy += y_v[j]
                sxy += x_v[j] * y_v[j]
                valid_n += 1
            end
        end
        if valid_n > 0
            mx = sx / valid_n
            my = sy / valid_n
            covs[n_neg + i] = (sxy / valid_n) - (mx * my)
        else
            covs[n_neg + i] = NaN
        end
    end
    return covs
end

function random_error_rannik(cx::Vector{Float64}, wy::Vector{Float64}, time_shift::Int64, sampling_span::Int64)
    #lagged_covs-> R_wc(Dong et al., 2021; Rannik, 2016) 

    lagged_covs = lagged_covariances(cx, wy, time_shift, sampling_span)
    
    random_error = stdf(lagged_covs)
    
    return random_error
end 
