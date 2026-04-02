#module OutputResults

#export output_results
#import ..Utilities as UT
#using ..Utilities
#using ..Config

#using DataFrames, CSV, Dates

function output_results(buffer)
    print("output_start...\n")
    ####set_colnames######
    colnames_timestamp = [:timestamp]
    colnames_wind = [:WD, :WS, :U, :V, :W, :W_sigma, :U_var, :V_var, :W_var, :UV_cov, :Ustar, :L]
    colnames_temperature = [:Ta, :T_sonic]
    colnames_gases = [:FC, :Rho_CO2, :Rho_H2O, :Rho_dry, :Rho_CO2_sigma, :Rho_H2O_sigma, :Mix_CO2, :Mix_H2O, :Mix_CO2_sigma, :Mix_H2O_sigma]
    colnames_heat = [:H, :LE]
      #colnames_lag = [:lag_FC, :lag_LE]
    colnames_qc = [:Stationarity, :random_error_fc, :n_co2_err, :n_csat_err]
    colnames_cp = [:T_cell, :p_cell]
    colnames_gases_2 = Symbol.(string.(colnames_gases) .* "_2")
    colnames_heat_2 = Symbol.(string.(colnames_heat) .* "_2")
      #colnames_lag_2 = Symbol.(string.(colnames_lag) .* "_2")
    colnames_qc_2 = Symbol.(string.(colnames_qc) .* "_2")
    colnames_cp_2 = Symbol.(string.(colnames_cp) .* "_2")
    ####timestamp#####
    df_timestamp = DataFrame([buffer.common.timestamp], colnames_timestamp)

    ####UNITS#########
    buffer.common.temperature.ta = buffer.common.temperature.ta .-273.15 #degC
    buffer.common.temperature.t_sonic =  buffer.common.temperature.t_sonic .-273.15 #degC
    buffer.irga1.gases.fc = buffer.irga1.gases.fc .* 1e6 #mu mol/m2/s
    buffer.irga1.gases.rho_co2 = buffer.irga1.gases.rho_co2 .* 1e3 #mmol/m3
    buffer.irga1.gases.rho_h2o = buffer.irga1.gases.rho_h2o .* 1e3 #mmol/m3
    buffer.irga1.gases.rho_co2_sigma = buffer.irga1.gases.rho_co2_sigma .* 1e3 #mmol/m3
    buffer.irga1.gases.rho_h2o_sigma = buffer.irga1.gases.rho_h2o_sigma .* 1e3 #mmol/m3
    buffer.irga1.gases.mix_co2 = buffer.irga1.gases.mix_co2 .* 1e6 #ppm
    buffer.irga1.gases.mix_h2o = buffer.irga1.gases.mix_h2o .* 1e3 #ppt
    buffer.irga1.gases.mix_co2_sigma = buffer.irga1.gases.mix_co2_sigma .* 1e6 #ppm
    buffer.irga1.gases.mix_h2o_sigma = buffer.irga1.gases.mix_h2o_sigma .* 1e3 #ppt
    lag_fc_1 = buffer.irga1.lag.lag_fc ./SAMPLING_SPAN #s
    lag_le_1 = buffer.irga1.lag.lag_le ./SAMPLING_SPAN #s
    #buffer.irga1.qc.stationarity = buffer.irga1.qc.stationarity * 100 #%
    buffer.irga1.qc.random_error_fc = buffer.irga1.qc.random_error_fc .* buffer.irga1.gases.rho_dry*1e6 #mu mol m-2 s-1
    if IS_IRGA1_CP
        buffer.cp1.t_cell = buffer.cp1.t_cell .-273.15 #degC
        buffer.cp1.p_cell = buffer.cp1.p_cell .* 1e-3 #kPa
        #print(buffer.cp1.t_cell[1])
    end

    if TWO_IRGA_FLAG
        buffer.irga2.gases.fc = buffer.irga2.gases.fc .* 1e6 #mu mol/m2/s
        buffer.irga2.gases.rho_co2 = buffer.irga2.gases.rho_co2 .* 1e3 #mmol/m3
        buffer.irga2.gases.rho_h2o = buffer.irga2.gases.rho_h2o .* 1e3 #mmol/m3
        buffer.irga2.gases.rho_co2_sigma = buffer.irga2.gases.rho_co2_sigma .* 1e3 #mmol/m3
        buffer.irga2.gases.rho_h2o_sigma = buffer.irga2.gases.rho_h2o_sigma .* 1e3 #mmol/m3
        buffer.irga2.gases.mix_co2 = buffer.irga2.gases.mix_co2 .* 1e6 #ppm
        buffer.irga2.gases.mix_h2o = buffer.irga2.gases.mix_h2o .* 1e3 #ppt
        buffer.irga2.gases.mix_co2_sigma = buffer.irga2.gases.mix_co2_sigma .* 1e6 #ppm
        buffer.irga2.gases.mix_h2o_sigma = buffer.irga2.gases.mix_h2o_sigma .* 1e3 #ppt
        lag_fc_2 = buffer.irga2.lag.lag_fc ./SAMPLING_SPAN #s
        lag_le_2 = buffer.irga2.lag.lag_le ./SAMPLING_SPAN #s
        #buffer.irga2.qc.stationarity = buffer.irga2.qc.stationarity * 100 #%
        buffer.irga2.qc.random_error_fc = buffer.irga2.qc.random_error_fc .* buffer.irga2.gases.rho_dry*1e6 #mu mol m-2 s-1
        if IS_IRGA2_CP
            buffer.cp2.t_cell = buffer.cp2.t_cell .-273.15 #degC
            buffer.cp2.p_cell = buffer.cp2.p_cell .* 1e-3 #kPa
        end
    end
    
    ####common########
    df_wind = df_from_struct(buffer.common.wind, colnames_wind)
    df_temperature = df_from_struct(buffer.common.temperature, colnames_temperature)   
    df_common = hcat(
        df_wind,
        df_temperature
    )
    
    ###IRGA1##########
    df_gases_irga1 = df_from_struct(buffer.irga1.gases, colnames_gases)
    df_heat_irga1 = df_from_struct(buffer.irga1.heat, colnames_heat)
    df_lag_irga1 = DataFrame(lag_fc = lag_fc_1, lag_le = lag_le_1)
    df_qc_irga1 = df_from_struct(buffer.irga1.qc, colnames_qc)
    
    df_irga1 = hcat(
        df_gases_irga1,
        df_lag_irga1,
        df_qc_irga1
    )

    if MONCRIEFF97_FLAG_1
        df_moncrieff_irga1 = DataFrame(factor_fc = buffer.irga1.moncrieff)
        df_irga1 = hcat(df_irga1, df_moncrieff_irga1)
    end
    if IS_IRGA1_CP
        df_cp_irga1 = df_from_struct(buffer.cp1, colnames_cp)
        df_irga1 = hcat(
            df_irga1,
            df_cp_irga1
        )
    end

    df_heat = df_heat_irga1
    ###IRGA2##########
    df_irga2 = nothing
    if TWO_IRGA_FLAG
        df_gases_irga2 = df_from_struct(buffer.irga2.gases, colnames_gases_2)
        df_heat_irga2 = df_from_struct(buffer.irga2.heat, colnames_heat_2)
        df_lag_irga2 = DataFrame(lag_fc_2 = lag_fc_2, lag_le_2 = lag_le_2)
        df_qc_irga2 = df_from_struct(buffer.irga2.qc, colnames_qc_2)
        
        df_irga2 = hcat(
        df_gases_irga2,
        df_lag_irga2,
        df_qc_irga2
        )
        
        if IS_IRGA2_CP
            df_cp_irga2 = df_from_struct(buffer.cp2, colnames_cp_2)
            df_irga2 = hcat(
            df_irga2,
            df_cp_irga2
            )
        end
        if MONCRIEFF97_FLAG_2
            df_moncrieff_irga2 = DataFrame(factor_fc_2=buffer.irga2.moncrieff)
            df_irga2 = hcat(df_irga2, df_moncrieff_irga2)
        end
        if ATMOSPHERIC_PARAMETER_IRGA == 1
            df_heat = DataFrame(H = df_heat_irga1.H, LE = df_heat_irga1.LE, LE_2 = df_heat_irga2.LE_2)
        else
            df_heat = DataFrame(H = df_heat_irga2.H_2, LE = df_heat_irga1.LE, LE_2 = df_heat_irga2.LE_2)
        end
    end
    
    ####summary#####
    df_summary = hcat(
        df_timestamp,
        df_common,
        df_irga1,
        df_heat
    )
    if TWO_IRGA_FLAG
        df_summary = hcat(
            df_summary,
            df_irga2
        )
    end
    filename = RESULT_OUT_DIRECTORY*"\\"*RESULT_FILE_NAME*".csv"
    CSV.write(filename, df_summary)
    print("output_end...\n")
    print("$(filename)\n")
end

#end