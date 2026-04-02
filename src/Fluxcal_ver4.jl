#using Statistics, Dates

#SRC_DIR = @__DIR__
#include(joinpath(SRC_DIR, "xcov_ver2.jl"))
#include(joinpath(SRC_DIR, "Obukhov_length_ver2.jl"))
#include(joinpath(SRC_DIR, "axisrotation_fix.jl"))
#include(joinpath(SRC_DIR, "winddirection_ver2.jl"))
#include(joinpath(SRC_DIR, "dsp_ver4.jl"))
#include(joinpath(SRC_DIR, "stationarity_test.jl"))
#include(joinpath(SRC_DIR, "rannik2016_ver2.jl"))

function compute_TTime_mean(timeseries)
    local df = DateFormat("yyyy-mm-dd HH:MM:SS:sss")
    local st = DateTime("0-12-31 0:00:00:0", df)
    return st + 
        Millisecond(
            round(
                meanf(
                    Dates.value.(
                        timeseries
                    )
                )
            )
        )
end

function compute_wind_base(u::Vector{Float64}, v::Vector{Float64}, w::Vector{Float64}; optrot::Int = 1)
    u_mean = meanf(u)
    v_mean = meanf(v)
    w_mean = meanf(w)
    wind_3d_mean = [u_mean, v_mean, w_mean]
    wd_offset = 0
    sonictype = "Campbell"
    wd = compute_wind_direction(wind_3d_mean, wd_offset, sonictype)
    #wind_3d_raw = hcat(u, v, w)
    ufx, ufy, ufz, theta, fai = axisrotation_modified(u, v, w, optrot)
    ws = meanf(sqrt.(ufx.^2 .+ ufy.^2))
    ustar = (covf(ufx, ufz)^2 + covf(ufy, ufz)^2)^0.25
    return (
        wd = wd,
        u_mean = u_mean,
        v_mean = v_mean,
        w_mean = w_mean,
        ws = ws,
        ustar = ustar,
        ufz = ufz
    )
end

function compute_rho_d(p::Union{Float64, Vector{Float64}}, t::Union{Vector{Float64}, Float64}, rho_h2o::Vector{Float64})
    return p ./ R ./ t .- rho_h2o
end

function compute_mixing_ratio(rho_x::Vector{Float64}, rho_d::Vector{Float64})
    return rho_x ./ rho_d
end

function compute_flux(
    coeff::Union{Int, Float64}, 
    w::Vector{Float64}, 
    x::Vector{Float64}; 
    lag_offset::Int64 = 0, 
    lag_min::Int64 = -20, 
    lag_max::Int64 = 10
)
    wx, lag = cross_cov(w, x, lag_offset, lag_min, lag_max)
    return coeff .* wx, lag
end

function compute_flux_wpl(
    rho_d_mean::Float64, 
    w::Vector{Float64}, 
    x::Vector{Float64}, 
    h2o::Vector{Float64}, 
    t_sonic::Vector{Float64}
)
    #ignoring the pressure term effect.
    #WPL correction
    rho_x_mean = meanf(x)
    rho_h2o_mean = meanf(h2o)
    rho_mean = rho_d_mean + rho_h2o_mean

    coeff_h2o = rho_x_mean / rho_d_mean
    coeff_t_sonic = rho_x_mean * (rho_mean / rho_d_mean) / meanf(t_sonic)
    
    wx, lag_offset_wpl = cross_cov(w, x, LAG_OFFSET, LAG_MIN, LAG_MAX)
    term_h2o = coeff_h2o * cross_cov(w, h2o, lag_offset_wpl, 0, 0)[1]
    term_t_sonic = coeff_t_sonic * covf(w, t_sonic)
    
    flux_x_op = wx + term_h2o + term_t_sonic
    
    return flux_x_op, lag_offset_wpl
end

function compute_fc_cp(
    rho_d_mean::Float64, 
    w::Vector{Float64}, 
    mix_co2::Vector{Float64}
)
    fc_cp, lag_fc = compute_flux(
        rho_d_mean, 
        w, 
        mix_co2, 
        lag_offset = LAG_OFFSET,
        lag_min = LAG_MIN,
        lag_max = LAG_MAX
    )
    return fc_cp, lag_fc
end

function compute_fc_op(
    rho_d_mean::Float64, 
    w::Vector{Float64}, 
    co2::Vector{Float64}, 
    h2o::Vector{Float64}, 
    t_sonic::Vector{Float64}
)
    fc_op, lag_fc = compute_flux_wpl(rho_d_mean, w, co2, h2o, t_sonic)
    return fc_op, lag_fc
end

function compute_le_cp(
    rho_d_mean::Float64, 
    w::Vector{Float64}, 
    mix_h2o::Vector{Float64}, 
    t_a_mean::Float64
)
    if t_a_mean >= T0
        lv = (2.50025 - 2365 * (t_a_mean - T0) * 1e-6) * 1e3
    else
        lv = (2.8341 - 149 * (t_a_mean - T0) * 1e-6) * 1e3
    end
    coeff_hl = rho_d_mean * lv * Mw
    le_cp, lag_le = compute_flux(
        coeff_hl, 
        w, 
        mix_h2o;
        lag_offset = LAG_OFFSET,
        lag_min = LAG_MIN,
        lag_max = LAG_MAX
    )
    return le_cp, lag_le
end

function compute_le_op(
    rho_d_mean::Float64, 
    w::Vector{Float64}, 
    h2o::Vector{Float64}, 
    t_sonic::Vector{Float64}, 
    t_a_mean::Float64
)
    if t_a_mean >= T0
        lv = (2.50025 - 2365 * (t_a_mean - T0) * 1e-6) * 1e3
    else
        lv = (2.8341 - 149 * (t_a_mean - T0) * 1e-6) * 1e3
    end
    lambda_wpl = lv * Mw
    e, lag_le = compute_flux_wpl(rho_d_mean, w, h2o, h2o, t_sonic)

    le_op = lambda_wpl *e
    return le_op, lag_le
end

function compute_specific_humidity(rho_h2o::Vector{Float64}, rho_d::Vector{Float64})
    num = rho_h2o * Mw
    den = rho_d .* Md .+ rho_h2o .* Mw
    return num ./ den #q: specific humidity (dimensionless)
end

function compute_h(
    rho_d::Vector{Float64}, 
    h2o::Vector{Float64}, 
    w::Vector{Float64}, 
    t_sonic::Vector{Float64}, 
    lag_le::Int
)
    #mole density
    rho_h2o_mean = meanf(h2o)
    rho_d_mean = meanf(rho_d)
    
    #moist air's heat capacity
    rho_cp = rho_d_mean * Md * Cpd + rho_h2o_mean * Mw * Cpv

    #specific humidity
    q = compute_specific_humidity(h2o, rho_d)
    q_mean = meanf(q)

    #alpha
    a = 0.51

    #first term: <w'Ts'>
    wts = covf(w, t_sonic)

    #second term:  <q> <w'Ts'>
    q_wts = q_mean * wts

    #third term: <Ts> <w'q'>
    ts_wq = meanf(t_sonic) * cross_cov(w, q, lag_le, 0, 0)[1]

    #summarized terms
    wta = wts - a * (q_wts +ts_wq)
    
    h = rho_cp * wta

    return h
end

function compute_t_a_mean(
    common_raw_calc::CommonRawData, 
    irga_raw_calc::IRGARawData, 
    p_atm::Float64 = PRESSURE_ATMOSPHERE
)
    rho_h2o_mean = meanf(irga_raw_calc.h2o)
    t_sonic_mean = meanf(common_raw_calc.t_sonic)
    rho_sonic = p_atm / (R * t_sonic_mean)
    mix_h2o_sonic = rho_h2o_mean / rho_sonic
    t_a_mean = t_sonic_mean / (1 + 0.32 * mix_h2o_sonic)
    return t_a_mean
end

function flux_calc_common(common_raw_calc::CommonRawData)
    u = dsp(common_raw_calc.u, 5, 10)
    v = dsp(common_raw_calc.v, 5, 10)
    w = dsp(common_raw_calc.w, 5, 10)
    t_sonic = common_raw_calc.t_sonic
    optrot = 1
    wind_base = compute_wind_base(u, v, w; optrot)
    wd = wind_base.wd
    u_mean = wind_base.u_mean
    v_mean = wind_base.v_mean
    w_mean = wind_base.w_mean
    ws = wind_base.ws
    ustar = wind_base.ustar
    w_rotated = wind_base.ufz
    t_sonic_mean = meanf(t_sonic) 
    wts = covf(w_rotated, t_sonic)
    L = obukhov(t_sonic_mean, PRESSURE_ATMOSPHERE, ustar, wts)
    w_sigma = stdf(w)
    u_var = varf(u)
    v_var = varf(v)
    w_var = varf(w)
    uv_cov = covf(u, v)

    wind = Wind(
        [wd], 
        [ws], 
        [u_mean], 
        [v_mean], 
        [w_mean],
        [w_sigma],
        [u_var], 
        [v_var],
        [w_var],
        [uv_cov],
        [ustar],
        [L]
    )
    
    return wind, [t_sonic_mean], w_rotated
end

function flux_calc_irga(
    common_raw_calc::CommonRawData, 
    irga_raw_calc::IRGARawData, 
    w_rotated::Vector{Float64}, 
    t_a_mean::Float64,
    is_cp_irga::Bool
)
    #@time begin
    #QC: NaN values
    n_csat_err = sum(isnan.(w_rotated))
    n_co2_err = sum(isnan.(irga_raw_calc.co2))
    
    #despike IRGA signals
    h2o = dsp(irga_raw_calc.h2o, 3.5, 10)
    co2 = dsp(irga_raw_calc.co2, 3.5, 10)

    #sonic temperature
    t_sonic = common_raw_calc.t_sonic

    #dry air density
    if is_cp_irga
        rho_d = compute_rho_d(irga_raw_calc.p_cell,irga_raw_calc.t_cell, h2o)
        #print(irga_raw_calc.p_cell[1])
    else
        rho_d = compute_rho_d(PRESSURE_ATMOSPHERE, t_a_mean, h2o)
    end
#end
    #@time begin
    #mixing ratios (vector)
    mix_co2 = compute_mixing_ratio(co2, rho_d)
    mix_h2o = compute_mixing_ratio(h2o, rho_d)
    #end
    #statistics
    #@time begin
    rho_d_mean = meanf(rho_d)
    rho_h2o_mean = meanf(h2o)
    rho_co2_mean = meanf(co2)
    rho_co2_sigma = stdf(co2)
    rho_h2o_sigma = stdf(h2o)
    mix_co2_mean = meanf(mix_co2)
    mix_h2o_mean = meanf(mix_h2o)
    mix_co2_sigma = stdf(mix_co2)
    mix_h2o_sigma = stdf(mix_h2o)
    #end
    #fluxes
    #@time begin
    if is_cp_irga
        fc, lag_fc = compute_fc_cp(rho_d_mean, w_rotated, mix_co2)
        le, lag_le = compute_le_cp(rho_d_mean, w_rotated, mix_h2o, t_a_mean)       
        stationarity_test = stationarity(mix_co2, w_rotated, lag_fc) 
    else
        fc, lag_fc = compute_fc_op(rho_d_mean, w_rotated, co2, h2o, t_sonic)
        le, lag_le = compute_le_op(rho_d_mean, w_rotated, h2o, t_sonic, t_a_mean)
        stationarity_test = stationarity(co2, w_rotated, lag_fc)
    end  
    #end   
    time_shift_max = 300 #sec
    random_error_fc = random_error_rannik(mix_co2, w_rotated, time_shift_max, SAMPLING_SPAN)
    h = compute_h(rho_d, h2o, w_rotated, t_sonic, lag_le)


    gases = Gases(
        [fc],
        [rho_co2_mean],
        [rho_h2o_mean],
        [rho_d_mean],
        [rho_co2_sigma],
        [rho_h2o_sigma],
        [mix_co2_mean],
        [mix_h2o_mean],
        [mix_co2_sigma],
        [mix_h2o_sigma]
    )

    heat = Heat(
        [h],
        [le]
    )

    lag = Lag(
        [lag_fc],
        [lag_le]
    )

    qc = QCResults(
        [stationarity_test],
        [random_error_fc],
        [n_co2_err],
        [n_csat_err]
    )

    return gases, heat, lag, qc
end

function cp(t_cell::Vector{Float64}, p_cell::Vector{Float64})
    return CPIRGAParameter([meanf(t_cell)], [meanf(p_cell)])
end