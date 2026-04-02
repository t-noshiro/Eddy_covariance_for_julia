#Author Taichi Noshiro
#using  Statistics, Interpolations


function mov(data::AbstractVector{Float64}, n::Int64, is_median::Bool)
    len = length(data)
    move = fill(NaN, len)
    #move = zeros(typeof(data[1]),length(data))
    #move = convert(Vector{Any}, move)
    half_high = div(n,2)
    half_low = rem(n, 2) == 1 ? half_high : half_high - 1

    tmp_buffer = Vector{Float64}(undef, n)
    for i = 1:len
        
        if i <= half_low
            s, e = 1, min(n, len)
        elseif i > len - half_high
            s, e = max(1, len - n + 1), len
        else
            s = max(1, i - half_low)
            e = min(len, i + half_high)
        end
        
        window_data = @view data[s:e]
        count = 0
        for val in window_data
            if !isnan(val)
                count += 1
                tmp_buffer[count] = val
            end
        end

        if count > 0
            valid_view = @view tmp_buffer[1:count]
            move[i] = is_median ? median(valid_view) : std(valid_view)
        end
        #move[i] = (opt == "omitnan") ? funcf(f, window_data, 1) : funcf(f, window_data, 0)
    end
    return move
end

movmedian(x, n) =  mov(x, n, true)
movstd(x, n) =  mov(x, n, false)

function  despike(x::AbstractVector, threshold_stds, n_window_width::Int, opt_dsp = 1)
    len_x = length(x)
    nseries = collect(1:len_x)
    xnew = copy(x)
    
    window_width = Int64(floor(len_x / n_window_width))
    mstdarray = movstd(x, window_width)
    mmedarray = movmedian(x, 3)

    xevaluate = zeros(len_x)
    fndnan = zeros(Bool, len_x)
    
    for i = 1:len_x
        i_end = min(i + 2, len_x)
        window_3 = @view x[i:i_end]
        mstd = mstdarray[i]
        mmed = mmedarray[i]
        if all(isnan.(window_3))
            fndnan[i] = 1
        else
            xevaluate[i] = threshold_stds * mstd - abs(x[i] - mmed)
        end
    end
    
    fndspikes = .!fndnan .& (xevaluate .< 0)
    num_spikes = sum(fndspikes)
    
    if num_spikes > length(x) -2
        println("caution:too many spikes!($num_spikes / $len_x)\n")
        return opt_dsp == 1 ? (x, num_spikes) : x
    end
    
    valid_mask = .!fndspikes .& .!fndnan
    
    if sum(valid_mask) < 2
        return opt_dsp == 1 ? (x, num_spikes) : x
    end
    
    itp_base = interpolate((nseries[valid_mask],), @view(x[valid_mask]), Gridded(Linear()))
    etp = extrapolate(itp_base, Line())

    if num_spikes > 0
        xnew[fndspikes] .= etp(nseries[fndspikes])
    end
    
    if opt_dsp == 1
        return xnew, num_spikes
    elseif opt_dsp == 2
        return xnew
    else
        return num_spikes
    end

end

dsp(x, threshold_stds, n_window_width) = despike(x, threshold_stds, n_window_width, 2)
