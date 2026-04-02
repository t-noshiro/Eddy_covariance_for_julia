#using Dates, DataFrames
mutable struct Gases
    fc ::Vector{Float64}
    rho_co2 ::Vector{Float64}
    rho_h2o ::Vector{Float64}
    rho_dry ::Vector{Float64}
    rho_co2_sigma ::Vector{Float64}
    rho_h2o_sigma ::Vector{Float64}
    mix_co2 ::Vector{Float64}
    mix_h2o ::Vector{Float64}
    mix_co2_sigma ::Vector{Float64}
    mix_h2o_sigma ::Vector{Float64}
end

mutable struct Wind
    wd ::Vector{Float64}
    ws ::Vector{Float64}
    u ::Vector{Float64}
    v ::Vector{Float64}
    w ::Vector{Float64}
    w_sigma ::Vector{Float64}
    u_var ::Vector{Float64}
    v_var ::Vector{Float64}
    w_var ::Vector{Float64}
    uv_cov ::Vector{Float64}
    ustar ::Vector{Float64}
    L ::Vector{Float64}
end

mutable struct Temperature
    ta ::Vector{Float64}
    t_sonic ::Vector{Float64}
end

mutable struct Heat
    h ::Vector{Float64}
    le ::Vector{Float64}
end

mutable struct Lag
    lag_fc ::Vector{Int64}
    lag_le ::Vector{Int64}
end

mutable struct QCResults
    stationarity ::Vector{Float64}
    random_error_fc ::Vector{Float64}
    n_co2_err ::Vector{Float64}
    n_csat_err ::Vector{Float64}
end

mutable struct FluxResultsCommon
    timestamp ::Vector{DateTime}
    wind ::Wind
    temperature ::Temperature
end

mutable struct FluxResultsIRGA
    timestamp ::Vector{DateTime}
    gases ::Gases
    heat ::Heat   
    lag ::Lag
    qc ::QCResults
    moncrieff ::Vector{Float64}
end

mutable struct CPIRGAParameter
    t_cell ::Vector{Float64}
    p_cell ::Vector{Float64}
end

mutable struct CommonRawData
    u ::Vector{Float64}
    v ::Vector{Float64}
    w ::Vector{Float64}
    t_sonic ::Vector{Float64}
end

mutable struct IRGARawData
    co2 ::Vector{Float64}
    h2o ::Vector{Float64}
    t_cell ::Union{Vector{Float64}, Nothing}
    p_cell ::Union{Vector{Float64}, Nothing}
end

mutable struct FluxBuffer
    common ::FluxResultsCommon
    irga1 ::FluxResultsIRGA
    irga2 ::FluxResultsIRGA
    cp1 ::Union{CPIRGAParameter, Nothing}
    cp2 ::Union{CPIRGAParameter, Nothing}
end

struct SeparationPeriod
    t_start ::DateTime
    t_end ::DateTime
    irga_separation ::Float64
end

function Gases_empty()
    return Gases(
        Float64[], #fc
        Float64[], #rho_co2
        Float64[], #rho_h2o
        Float64[], #rho_dry
        Float64[], #rho_co2_sigma
        Float64[], #rho_h2o_sigma
        Float64[], #mix_co2
        Float64[], #mix_h2o
        Float64[], #mix_co2_sigma
        Float64[] #mix_h2o_sigma
    )
end

function Wind_empty()
    return Wind(
        Float64[],  # wd
        Float64[],  # ws
        Float64[],  # u
        Float64[],  # v
        Float64[],  # w
        Float64[],  # w_sigma
        Float64[],  # u_var
        Float64[],  # v_var
        Float64[],  # w_var
        Float64[],  # uv_cov
        Float64[],  # ustar
        Float64[]   # L
    )
end

function Temperature_empty()
    return Temperature(
        Float64[],  # ta
        Float64[]  # t_sonic
    )
end

function Heat_empty()
    return Heat(
        Float64[],
        Float64[]
    )
end

function Lag_empty()
    return Lag(
        Int[],  # lag_fc
        Int[]   # lag_le
    )
end

function QCResults_empty()
    return QCResults(
        Float64[],  # stationarity
        Float64[],  # random_error_fc
        Float64[],  # n_co2_err
        Float64[]   # n_csat_err
    )
end

function FluxResultsCommon_empty()
    return FluxResultsCommon(
        DateTime[],      # timestamp
        Wind_empty(),    # wind
        Temperature_empty()   # heat
    )
end

function FluxResultsIRGA_empty()
    return FluxResultsIRGA(
        DateTime[],      # timestamp
        Gases_empty(), #gases
        Heat_empty(),   # le
        Lag_empty(), #lag
        QCResults_empty(), #qc
        Float64[]   # moncrieff
    )
end

function CPIRGAParameter_empty()
    return CPIRGAParameter(
        Float64[],
        Float64[]
    )
end

function CommonRawData_empty()
    return CommonRawData(
        Float64[],
        Float64[],
        Float64[],
        Float64[]
    )
end

function IRGARawData_empty(cp_flag)
    if cp_flag
        return IRGARawData(
        Float64[],
        Float64[],
        Float64[],
        Float64[]
        )
    else
        return IRGARawData(
            Float64[],
            Float64[],
            nothing,
            nothing
        )
    end
end

function FluxBuffer_empty()
    return FluxBuffer(
        FluxResultsCommon_empty(),
        FluxResultsIRGA_empty(),
        FluxResultsIRGA_empty(),
        IS_IRGA1_CP ? CPIRGAParameter_empty() : nothing,
        IS_IRGA2_CP ? CPIRGAParameter_empty() : nothing
    )
end

function append_common_raw!(raw::CommonRawData, tmp::CommonRawData)
    append!(raw.u, tmp.u)
    append!(raw.v, tmp.v)
    append!(raw.w, tmp.w)
    append!(raw.t_sonic, tmp.t_sonic)
end

function append_irga_raw!(raw::IRGARawData, tmp::IRGARawData)
    append!(raw.co2, tmp.co2)
    append!(raw.h2o, tmp.h2o)
    if !isnothing(raw.t_cell)
        append!(raw.t_cell, tmp.t_cell)
        append!(raw.p_cell, tmp.p_cell)
    end
end

function append_wind!(res ::Wind, calc_res ::Wind)
    append!(res.wd, calc_res.wd)
    append!(res.ws, calc_res.ws)
    append!(res.u, calc_res.u)
    append!(res.v, calc_res.v)
    append!(res.w, calc_res.w)
    append!(res.w_sigma, calc_res.w_sigma)
    append!(res.u_var, calc_res.u_var)
    append!(res.v_var, calc_res.v_var)
    append!(res.w_var, calc_res.w_var)
    append!(res.uv_cov, calc_res.uv_cov)
    append!(res.ustar, calc_res.ustar)
    append!(res.L, calc_res.L)
end

function append_temperature!(res ::Temperature, calc_res::Temperature)
    append!(res.ta, calc_res.ta)
    append!(res.t_sonic, calc_res.t_sonic)
end

function append_heat!(res::Heat, calc_res::Heat)
    append!(res.h, calc_res.h)
    append!(res.le, calc_res.le)
end

function append_gases!(res ::Gases, calc_res ::Gases)
    append!(res.fc, calc_res.fc)
    append!(res.rho_co2, calc_res.rho_co2)
    append!(res.rho_h2o, calc_res.rho_h2o)
    append!(res.rho_dry, calc_res.rho_dry)
    append!(res.rho_co2_sigma, calc_res.rho_co2_sigma)
    append!(res.rho_h2o_sigma, calc_res.rho_h2o_sigma)
    append!(res.mix_co2, calc_res.mix_co2)
    append!(res.mix_h2o, calc_res.mix_h2o)
    append!(res.mix_co2_sigma, calc_res.mix_co2_sigma)
    append!(res.mix_h2o_sigma, calc_res.mix_h2o_sigma)
end

function append_lag!(res ::Lag, calc_res::Lag)
    append!(res.lag_fc, calc_res.lag_fc)
    append!(res.lag_le, calc_res.lag_le)
end

function append_qc!(res ::QCResults, calc_res ::QCResults)
    append!(res.n_co2_err, calc_res.n_co2_err)
    append!(res.n_csat_err, calc_res.n_csat_err)
    append!(res.stationarity, calc_res.stationarity)
    append!(res.random_error_fc, calc_res.random_error_fc)
end

function append_cp_irga_parameter!(res ::CPIRGAParameter, calc_res ::CPIRGAParameter)
    append!(res.t_cell, calc_res.t_cell)
    append!(res.p_cell, calc_res.p_cell)
end

function divide_vector(target_vec, time_span, sampling_span)
    segment_end_index = Int64(time_span * 60 * sampling_span)

    if length(target_vec) > segment_end_index
        remaining_segment = target_vec[segment_end_index + 1:end]
        calc_segment = target_vec[1:segment_end_index]

    elseif length(target_vec) == segment_end_index
        remaining_segment = zeros(typeof(target_vec[1]), 0)
        calc_segment = target_vec[1:segment_end_index]

    else
        remaining_segment = target_vec
        calc_segment = nothing
    end
    
    return calc_segment, remaining_segment
end

function divide_common(common_raw ::CommonRawData, time_span, sampling_span)
    u_calc, u_rem = divide_vector(common_raw.u, time_span, sampling_span)
    v_calc, v_rem = divide_vector(common_raw.v, time_span, sampling_span)
    w_calc, w_rem = divide_vector(common_raw.w, time_span, sampling_span)
    t_sonic_calc, t_sonic_rem = divide_vector(common_raw.t_sonic, time_span, sampling_span)
    if isnothing(u_calc)
        return CommonRawData_empty(), CommonRawData(u_rem, v_rem, w_rem, t_sonic_rem)
    else
        return CommonRawData(u_calc, v_calc, w_calc, t_sonic_calc), 
        CommonRawData(u_rem, v_rem, w_rem, t_sonic_rem)
    end
end

function divide_irga(irga_raw ::IRGARawData, time_span, sampling_span)
    co2_calc, co2_rem = divide_vector(irga_raw.co2, time_span, sampling_span)
    h2o_calc, h2o_rem = divide_vector(irga_raw.h2o, time_span, sampling_span)
    if !isnothing(irga_raw.t_cell)
        t_cell_calc, t_cell_rem = divide_vector(irga_raw.t_cell, time_span, sampling_span)
        p_cell_calc, p_cell_rem = divide_vector(irga_raw.p_cell, time_span, sampling_span)
    else
        t_cell_calc, t_cell_rem = nothing, nothing
        p_cell_calc, p_cell_rem = nothing, nothing
    end

    #This part is following the form, "return calc_segment, remaining_segment"
    if isnothing(co2_calc)
        return IRGARawData_empty(false), irga_raw #This calc_segment is gabage. DON'T USE THIS CALC_SEGMENT.
    else
        return IRGARawData(co2_calc, h2o_calc, t_cell_calc, p_cell_calc),
        IRGARawData(co2_rem, h2o_rem, t_cell_rem, p_cell_rem)
    end
end

function df_from_struct(obj, colnames::Vector{Symbol})
    fields = fieldnames(typeof(obj))
    @assert length(fields) == length(colnames) "Number of fields and column names do not match"

    cols = [getfield(obj, f) for f in fields]
    return DataFrame(colnames .=> cols)
end

function irga_separation_for(TTime::DateTime, table)
    for p in table
        if p.t_start <= TTime < p.t_end
            return p.irga_separation
        end
    end
end