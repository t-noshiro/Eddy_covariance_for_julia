# Stationarity for covariance and "dispersion (variance, i=j)"
# based on Foken & Wichura (1996, AFM)
# FW96 suggests Stationarity to be less than 30% (P94).
# Author Hiroki Ikawa
# Edited by Taichi Noshiro
function stationarity(xin::Vector{Float64},yin::Vector{Float64},lag::Int64)
    ndiv = 6 # the number of segments
    #lag = Int64(lag)
    len_in = length(xin)

    s_x, e_x = max(1, 1 + lag), min(len_in, len_in + lag)
    s_y, e_y = max(1, 1 - lag), min(len_in, len_in - lag)
    @views begin
        x = xin[s_x : e_x]
        y = yin[s_y : e_y]
        ndata = length(x)
        nseg = floor(Int, ndata / ndiv)
     
        covlegt = 0.0
        valid_divs = 0 
        for i = 1:ndiv
            idx_start = nseg * (i - 1) + 1
            idx_end = nseg * i
            xleg = x[idx_start : idx_end]
            yleg = y[idx_start : idx_end]
            covleg = covf(xleg, yleg) # Eq.13
            covlegt = covlegt + covleg
            if !isnan(covleg)
                valid_divs += 1
            end
        end 
     
        if valid_divs == 0
            return NaN
        end
        covlegt = covlegt / valid_divs  # Eq.14
        covrun = covf(x, y)       # Eq.15
 
        stationarity = 100. * abs((covrun - covlegt)/covrun)
    end 
    return stationarity
end