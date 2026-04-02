## This program returns wind direction, wd based on wind_3d_mean and sonic type.
# wind_3d_mean[1]: mean Ux, wind_3d_mean[2]: mean Uy, wind_3d_mean[3]: mean Uz
# sonic type: "Gill", "Campbell"
#Author Hiroki Ikawa
#Edited by Docopy

function compute_wind_direction(wind_3d_mean ::Vector{Float64}, wd_offset, sonictype ::String)
    u_mean = wind_3d_mean[1]
    v_mean = wind_3d_mean[2]
    wind_angle_rad = atan(abs(v_mean)/abs(u_mean))
    wind_angle = wind_angle_rad / 2 / pi * 360 #deg

    if sonictype == "Gill"
    
        if (u_mean >0 && v_mean > 0) 
            wd = 180 - wind_angle
        elseif (u_mean > 0 && v_mean < 0) 
            wd = 180 + wind_angle
        elseif (u_mean < 0 && v_mean > 0) 
            wd = wind_angle
        elseif (u_mean < 0 && v_mean < 0) 
            wd = 360 - wind_angle
        end

    elseif sonictype == "Campbell"

        if (u_mean > 0 && v_mean > 0)
            wd = 360 - wind_angle
        elseif (u_mean > 0 && v_mean < 0)
            wd = wind_angle
        elseif (u_mean < 0 && v_mean > 0)
            wd = 180 + wind_angle
        elseif (u_mean < 0 && v_mean < 0) 
            wd = 180 - wind_angle
        else
            wd = NaN
        end    
    
    end

    wd = wd + wd_offset

    if wd < 0
        wd = wd + 360 
    end

    if wd > 360
        wd = wd - 360
    end


    return wd

end