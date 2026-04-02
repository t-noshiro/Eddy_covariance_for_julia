# Obkhov length calculation following P44 of Foken (2008)
# Author Hiroki Ikawa
# Edited by Taichi Noshiro
# 20260108 changed the units to SI standards and variable names.
# Note that for EC calculation, I use sonic temperature as Tv and 
# covariance of sonic temp and vertical wind for wTv, though 
# virtual temp and sonic temp are slightly different P144 (Foken)
# Input (Unit is following the SI standard.)
# t_v      virtual temperature (K) 
# p_atm   atmospheric pressure (Pa)
# ustar   friction velocity (m/s)
# wtv     flux of virtual temperature (K*m/s)

# Output 
# L       Obkhov length (m)

function obukhov(t_v_mean ::Float64, p_atm::Float64, ustar::Float64, wTv::Float64)
    t_pot = t_v_mean* ((1e5 / p_atm)^0.286) # Potential temperature (K)
    L = -(ustar) ^ 3 / (0.4 * 9.8 / t_pot * wTv)
    return L
end