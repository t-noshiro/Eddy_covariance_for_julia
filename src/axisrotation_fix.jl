###########################################################################
#Author Taichi Noshiro
## INPUT 
# u, v, w: measured wind speed of x-, y- and z-axis components from sonic anemometer, respectively (m s-1)
## OUTPUT
# u_f2, v_f2, w_f2: rotated wind speed (m s-1)
# If you want to apply only the first rotation, then change the optrot from 2 to 1
###########################################################################
function axisrotation_modified(u::AbstractVector{Float64}, v::AbstractVector{Float64}, w::AbstractVector{Float64}, optrot::Int)

    # 2-D rotations 
    u_mean = meanf(u)
    v_mean = meanf(v)
    theta = atan(v_mean / u_mean)
    fai = 0.0
    
    rot1 = [
        cos(theta) sin(theta) 0
        -sin(theta) cos(theta) 0
        0 0 1
    ]

    u1 = rot1 * transpose([u v w])
    u1 = transpose(u1)
    #=
    c1, s1 = cos(theta), sin(theta)
    u_f1 = @. u * c1 + v * s1
    v_f1 = @. -s1 * u + c1 * v
    w_f1 = w
    =#
    if optrot == 1
        #return u_f1, v_f1, w_f1, theta, fai
        return u1[:, 1], u1[:, 2], u1[:,3], theta, fai
    elseif optrot == 2

        fai = atan(meanf(w_f1) / meanf(u_f1))
        #c2, s2 = cos(fai), sin(fai)
        
        rot2 = [
            cos(fai) 0 sin(fai)
            0 1 0
            -sin(fai) 0 cos(fai)
        ]

        u2 = rot2 * transpose(u1)
        u2 = transpose(u2)
        

        #u_f2 = @. c2 * u_f1 + s2 * w_f1
        #v_f2 = @. v_f1
        #w_f2 = @. -s2 * u_f1 + c2 * w_f1

        #return u_f2, v_f2, w_f2, theta, fai
        return u2[:,1], u2[:,2], u2[:,3], theta, fai
    end

    

end