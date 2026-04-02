#########################
#This version is NOT usable for old ver ECsystems. (This message is for "Author".)
#Author: Docopypl (git)
#########################
#module MakingData_v5

#export search_file, make_raw

#import ..Utilities as UT
#using ..Config

#using CSV
#using DataFrames
#using Dates

function search_file(data_dir::String, keyword::String)
    files = readdir(data_dir; join = true)
    target = filter(f->occursin(keyword, lowercase(f)), files)

    if isempty(target)
        error("Now raw data file found in $(data_dir) matching keyword: $(keyword)")
    end

    return target
end

function make_common_raw(data)
    u = data[:, U_COLUMN_NAME]
    v = data[:, V_COLUMN_NAME]
    w = data[:, W_COLUMN_NAME]
    if @isdefined SOS_COLUMN_is
        if SOS_COLUMN_is == "Ts"
            t_sonic = data[:, SOS_COLUMN_NAME] .+ T0 #K
        elseif SOS_COLUMN_is == "sos"
            sos = data[:, SOS_COLUMN_NAME]
            sos = dsp(sos, 5, 10)
            t_sonic = (sos.^2)./1.4./287.04
        end
    else
        t_sonic = data[:,SOS_COLUMN_NAME]
    end
    return CommonRawData(u, v, w, t_sonic)
end

function make_irga1_raw(data)
    co2 = data[:, CO2_DENSITY_COLUMN_NAME] .* CONVERT_UNIT_CO2_DENSITY
    h2o = data[:, WATER_VAPOUR_DENSITY_COLUMN_NAME] .* CONVERT_UNIT_WATER_VAPOUR_DENSITY
    t_cell = nothing
    p_cell = nothing
    if IS_IRGA1_CP
        t_in = data[:, T_IN_COLUMN_NAME]
        t_out = data[:, T_OUT_COLUMN_NAME]
        t_cell = 0.8.*t_in .+ 0.2.*t_out .+ T0
        p_cell = data[:,PRESSURE_COLUMN_NAME] .* CONVERT_UNIT_PRESSURE
    end
    return IRGARawData(co2, h2o, t_cell, p_cell)
end

function make_irga2_raw(data)
    co2 = data[:, ANOTHER_CO2_DENSITY_COLUMN_NAME] .* CONVERT_UNIT_ANOTHER_CO2_DENSITY
    h2o = data[:, ANOTHER_WATER_VAPOUR_DENSITY_COLUMN_NAME] .* CONVERT_UNIT_ANOTHER_WATER_VAPOUR_DENSITY
    t_cell = nothing
    p_cell = nothing
    if IS_IRGA2_CP
        t_in = data[:, ANOTHER_T_IN_COLUMN_NAME]
        t_out = data[:, ANOTHER_T_OUT_COLUMN_NAME]
        t_cell = 0.8.*t_in .+ 0.2.*t_out .+ T0
        p_cell = data[:,ANOTHER_PRESSURE_COLUMN_NAME] .* CONVERT_UNIT_ANOTHER_PRESSURE
    end
    return IRGARawData(co2, h2o, t_cell, p_cell)
end

function make_raw(file::String, TWO_IRGA_FLAG::Bool)
    select_cols = [
        TIME_COLUMN_NAME, 
        CO2_DENSITY_COLUMN_NAME, 
        WATER_VAPOUR_DENSITY_COLUMN_NAME, 
        U_COLUMN_NAME, 
        V_COLUMN_NAME, 
        W_COLUMN_NAME, 
        SOS_COLUMN_NAME, 
        PRESSURE_COLUMN_NAME, 
        T_IN_COLUMN_NAME, 
        T_OUT_COLUMN_NAME
    ]
    if TWO_IRGA_FLAG
        append!(
            select_cols, 
            [
                ANOTHER_CO2_DENSITY_COLUMN_NAME, 
                ANOTHER_WATER_VAPOUR_DENSITY_COLUMN_NAME, 
                ANOTHER_PRESSURE_COLUMN_NAME, 
                ANOTHER_T_IN_COLUMN_NAME, 
                ANOTHER_T_OUT_COLUMN_NAME
            ]
        )
    end
    select_cols = Symbol.(select_cols)
    data = CSV.read(
        file, 
        DataFrame; 
        header = HEADER_LINE, 
        skipto = START_LINE, 
        delim=",", 
        select=select_cols
    )
    date_format = DateFormat(DATE_FORMAT*SEPARATING_CHARACTER*TIME_FORMAT)
    if IS_SEPARATING_DATE_AND_TIME
        Day = Dates.format.(data[:,DATE_COLUMN_NAME], DATE_FORMAT)
        TIME_tmp = data[:, TIME_COLUMN_NAME]
        TTime_str = Day.*SEPARATING_CHARACTER .* TIME_tmp
    else
        TTime_str = data[:,TIME_COLUMN_NAME]
    end

    TTime = DateTime.(TTime_str, date_format)

    buf_common_raw = make_common_raw(data)
    buf_irga1_raw = make_irga1_raw(data)
    if TWO_IRGA_FLAG
        buf_irga2_raw = make_irga2_raw(data)
        return TTime, buf_common_raw, buf_irga1_raw, buf_irga2_raw
    end
    return TTime, buf_common_raw, buf_irga1_raw
end