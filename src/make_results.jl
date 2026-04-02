
#module MakingResults

#export make_results

#using ..Config
#SRC_DIR = @__DIR__

#import ..Utilities as UT
#using ..Utilities


#include(joinpath(SRC_DIR,"Fluxcal_ver4.jl"))
#include(joinpath(SRC_DIR,"makingdata_ver5.0.jl"))
#include(joinpath(SRC_DIR,"Moncrieff97.jl"))
#using Dates
#using ..MakingData_v5

function process_segment(
    buffer, timestamp, common_raw_calc, irga1_raw_calc, irga2_raw_calc,
    temperature, table1, table2
)
    wind,
    temperature.t_sonic,
    w_rotated = flux_calc_common(common_raw_calc)

    if TWO_IRGA_FLAG
        if !IS_IRGA1_CP
            temperature.ta = [compute_t_a_mean(common_raw_calc, irga1_raw_calc)]
        elseif !IS_IRGA2_CP
            temperature.ta = [compute_t_a_mean(common_raw_calc, irga2_raw_calc)]
        else
            temperature.ta = [temperature.t_sonic]
        end
    else
        if !IS_IRGA1_CP
            temperature.ta = [compute_t_a_mean(common_raw_calc, irga1_raw_calc)]
        else
            temperature.ta = [temperature.t_sonic]
        end
    end
        
    print("irga1_flux\n")
    gases_irga1,
    heat_irga1,
    lag_irga1,
    qc_results_irga1 = flux_calc_irga(
        common_raw_calc,
        irga1_raw_calc,
        w_rotated,
        temperature.ta[1],
        IS_IRGA1_CP
    )
    print("irga2_flux\n")
    if TWO_IRGA_FLAG
        gases_irga2,
        heat_irga2,
        lag_irga2,
        qc_results_irga2 = flux_calc_irga(
        common_raw_calc,
        irga2_raw_calc,
        w_rotated,
        temperature.ta[1],
        IS_IRGA2_CP
        )
    end

    if IS_IRGA1_CP
        cp1_tmp = cp(irga1_raw_calc.t_cell, irga1_raw_calc.p_cell)
    end

    if IS_IRGA2_CP
        cp2_tmp = cp(irga2_raw_calc.t_cell, irga2_raw_calc.p_cell)
    end

    print("moncrieff\n") 
    if MONCRIEFF97_FLAG_1 && !isnothing(table1)  
        separation_irga1 = irga_separation_for(timestamp, table1)
        _, _, _, factor_fc = Moncrieff97(
            Z_OBSERVE_1,
            wind.ws[1],
            wind.L[1],
            separation_irga1,
            TUBE_RADIUS_1,
            TUBE_LENGTH_1,
            PUMP_FLOW_1,
            Float64(SAMPLING_SPAN),
            SAT_PATH_LEN,
            IRGA_PATH_LEN_1
        )
        factor_fc = factor_fc * ADDITIONAL_FACTOR
        gases_irga1.fc = gases_irga1.fc .* factor_fc
        moncrieff_irga1 = factor_fc
    end

     
    if TWO_IRGA_FLAG && MONCRIEFF97_FLAG_2   
        separation_irga2 = irga_separation_for(timestamp, table2)
        _, _, _, factor_fc = Moncrieff97(
            Z_OBSERVE_2,
            wind.ws[1],
            wind.L[1],
            separation_irga2,
            TUBE_RADIUS_2,
            TUBE_LENGTH_2,
            PUMP_FLOW_2,
            Float64(SAMPLING_SPAN),
            SAT_PATH_LEN,
            IRGA_PATH_LEN_2
        )
        factor_fc = factor_fc * ADDITIONAL_FACTOR
        gases_irga2.fc = gases_irga2.fc .* factor_fc
        moncrieff_irga2 = factor_fc
    end


    print(timestamp, "\n")
    #############
    #data-output#
    #############
    push!(buffer.common.timestamp, timestamp)
    append_wind!(buffer.common.wind, wind)
    append_temperature!(buffer.common.temperature, temperature)
    append_gases!(buffer.irga1.gases, gases_irga1)
    append_heat!(buffer.irga1.heat, heat_irga1)
    append_lag!(buffer.irga1.lag, lag_irga1)
    append_qc!(buffer.irga1.qc, qc_results_irga1)

    if MONCRIEFF97_FLAG_1
        push!(buffer.irga1.moncrieff, moncrieff_irga1)
    end

    if TWO_IRGA_FLAG
        append_gases!(buffer.irga2.gases, gases_irga2)
        append_heat!(buffer.irga2.heat, heat_irga2)
        append_lag!(buffer.irga2.lag, lag_irga2)
        append_qc!(buffer.irga2.qc, qc_results_irga2)

        if MONCRIEFF97_FLAG_2
            push!(buffer.irga2.moncrieff, moncrieff_irga2)
        end
        print(cp1_tmp.t_cell[1])
        if IS_IRGA1_CP
            append_cp_irga_parameter!(buffer.cp1, cp1_tmp)
        end

        if IS_IRGA2_CP
            append_cp_irga_parameter!(buffer.cp2, cp2_tmp)
        end
    end
end


function make_results(file_list, buffer)
    n_files = length(file_list)
    
    ############
    #initialize#
    ############
    TTime = DateTime[] 
    common_raw = CommonRawData_empty()
    irga1_raw = IRGARawData_empty(IS_IRGA1_CP)
    irga2_raw = TWO_IRGA_FLAG ? IRGARawData_empty(IS_IRGA2_CP) : nothing
    IRGA_SEPARATION_TABLE_1 = nothing
    IRGA_SEPARATION_TABLE_2 = nothing

    if MONCRIEFF97_FLAG_1
        dataformat_sep = DateFormat("yyyy-mm-ddTHH:MM:SS")
        IRGA_SEPARATION_TABLE_1 = Vector{Union{SeparationPeriod, Nothing}}(nothing, length(IRGA_SEPARATIONS_1))
        for i in 1:length(IRGA_SEPARATIONS_1)
            IRGA_SEPARATION_TABLE_1[i] = SeparationPeriod(
                DateTime(PERIOD_STARTS[i], dataformat_sep), 
                DateTime(PERIOD_ENDS[i], dataformat_sep), 
                IRGA_SEPARATIONS_1[i]
            )
        end
        if TWO_IRGA_FLAG && MONCRIEFF97_FLAG_2
            IRGA_SEPARATION_TABLE_2 = Vector{Union{SeparationPeriod, Nothing}}(nothing, length(IRGA_SEPARATIONS_2))
            for i in 1:length(IRGA_SEPARATIONS_2)
                IRGA_SEPARATION_TABLE_2[i] = SeparationPeriod(
                    DateTime(PERIOD_STARTS[i], dataformat_sep),
                    DateTime(PERIOD_ENDS[i], dataformat_sep), 
                    IRGA_SEPARATIONS_2[i]
                )
            end
        end
    end

    i = 1
  
    while i <= n_files
        ############
        #makingdata#
        ############
        data_file = file_list[i]

        if TWO_IRGA_FLAG
            TTime_tmp,
            common_raw_tmp, 
            irga1_raw_tmp, 
            irga2_raw_tmp = make_raw(data_file, TWO_IRGA_FLAG)
        else
            TTime_tmp, 
            common_raw_tmp,
            irga1_raw_tmp = make_raw(data_file, TWO_IRGA_FLAG)
        end
      
        print("################################\n")
        print(i, "/$(n_files)_read_compl\n")
        #print(length(TTime_local), "\n")

        ################
        #attaching_data#
        ################
        if length(TTime) == 0||TTime[end] - TTime_tmp[1] > Millisecond(floor(1/SAMPLING_SPAN*1000) + 1)
            TTime = copy(TTime_tmp)
            common_raw = deepcopy(common_raw_tmp)
            irga1_raw = deepcopy(irga1_raw_tmp)
            if TWO_IRGA_FLAG
                irga2_raw = deepcopy(irga2_raw_tmp)
            end
        else
            append!(TTime, TTime_tmp)
            append_common_raw!(common_raw, common_raw_tmp)
            append_irga_raw!(irga1_raw, irga1_raw_tmp)
            if TWO_IRGA_FLAG
                append_irga_raw!(irga2_raw, irga2_raw_tmp)
            end
        end
    
        ##########################
        #cutting_data&calculating#
        ##########################
        #@time begin
        print("cutting\n")
        while true
            #@time begin
            TTime_calc, TTime_rem = divide_vector(TTime, TIME_SPAN, SAMPLING_SPAN)
            common_raw_calc, common_raw_rem = divide_common(common_raw, TIME_SPAN, SAMPLING_SPAN)
            irga1_raw_calc, irga1_raw_rem = divide_irga(irga1_raw, TIME_SPAN, SAMPLING_SPAN)
            if TWO_IRGA_FLAG
                irga2_raw_calc, irga2_raw_rem = divide_irga(irga2_raw, TIME_SPAN, SAMPLING_SPAN)
            end

            if isnothing(TTime_calc)
                break
            end
            #end
            ##################
            #Flux_calculation#
            ##################
            #@time begin
            print("cal_executing...\n")
            print("common\n")
            temperature = Temperature_empty()
        
            timestamp = compute_TTime_mean(TTime_calc)
            #end    
            process_segment(buffer, timestamp, common_raw_calc, irga1_raw_calc, irga2_raw_calc,
                temperature, IRGA_SEPARATION_TABLE_1, IRGA_SEPARATION_TABLE_2
            ) 

            TTime = TTime_rem
            common_raw = common_raw_rem
            irga1_raw = irga1_raw_rem
            if TWO_IRGA_FLAG
                irga2_raw = irga2_raw_rem
            end
            print("calc_compl...\n")
         #end
        end
        i += 1
    end
end

#end