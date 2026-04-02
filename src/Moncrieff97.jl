#Author Hiroki Ikawa
# High frequency correction (Moncrieff97, Moore86) in Julia

#=
# INPUT (example values)
zobs = 2.0                   # measurement height [m]
U = 1.0                      # wind speed [m/s]
L = -5.0                     # Obukhov length [m]

IRGA_SEPARATION = 0.3        # sensor separation [m]
mHz = 10.0                   # sampling frequency [Hz]
TUBE_RADIUS = 2e-3           # tube radius [m]
TUBE_LENGTH = 0.5            # tube length [m]
PUMP_FLOW = 1.5e-4           # flow rate [m3 s-1]

SAT_PATH_LEN = 0.1                   # sonic path length [m]
IRGA_PATH_LEN = 0.13               # LI-COR path length [m]
=#

function Moncrieff97(zobs::Float64, U::Float64, L::Float64,
    IRGA_SEPARATION::Float64,
    TUBE_RADIUS::Float64, TUBE_LENGTH::Float64, PUMP_FLOW::Float64,   
    mHz::Float64,
    SAT_PATH_LEN::Float64 = 0.1, # value for CSAT3
    IRGA_PATH_LEN::Float64 = 0.13 # value for LI7200 and LI7500
    )


z = zobs
tauq = 1.0 / mHz   # 
LF = 2.073
C = 4.115

# Molecular diffusion coefficient for water vapor and CO2
Dq = 2.4e-5
Dc = 1.4e-5

n = zeros(Float64, 18)
n[1] = 1e-5

Twu = ones(Float64, 18)
Tww = ones(Float64, 18)
Twq = ones(Float64, 18)
Twc = ones(Float64, 18)
Twt = ones(Float64, 18)
Ts  = ones(Float64, 18)
Tdq  = ones(Float64, 18)
Tdc  = ones(Float64, 18)
Ttq  = ones(Float64, 18)
Ttc  = ones(Float64, 18)

UW = zeros(Float64, 18)
WX = zeros(Float64, 18)
simp = zeros(Float64, 18)

Suw = zeros(Float64, 18)
Swt = zeros(Float64, 18)
Swq = zeros(Float64, 18)
Swc = zeros(Float64, 18)

for i in 1:18
    # Line averaging
    flineSonic = n[i] * SAT_PATH_LEN / U
    xSonic = 2π * flineSonic

    flineIrga = n[i] * IRGA_PATH_LEN / U
    xIrga = 2π * flineIrga

    if xSonic > 0.02
        Twu[i] = (1.0 / xSonic) * (3 + exp(-xSonic) - 4 * (1 - exp(-xSonic)) / xSonic)
    end

    if xIrga > 0.02
        tmp = (1.0 / xIrga) * (3 + exp(-xIrga) - 4 * (1 - exp(-xIrga)) / xIrga)
        Twc[i] = tmp
        Twq[i] = tmp
    end

    if xSonic > 0.04
        Tww[i] = (4.0 / xSonic) * (1 + exp(-xSonic)/2 - 3 * (1 - exp(-xSonic)) / (2 * xSonic))
    end

    Twt[i] = Twu[i]

    # Sensor separation
    fsep = n[i] * IRGA_SEPARATION / U
    if fsep > 0.01
        Ts[i] = exp(-9.9 * fsep^1.5)
    end

    # Time constant gain
    if tauq > 0
        Tdq[i] = 1 ./ sqrt(1 + (2π * n[i] * tauq)^2)
        Tdc[i] = 1 ./ sqrt(1 + (2π * n[i] * tauq)^2)
    end

    Ttq[i] = exp(-pi^2*TUBE_RADIUS^2*n[i]^2*pi*TUBE_RADIUS^2*TUBE_LENGTH/PUMP_FLOW/6/Dq)
    Ttc[i] = exp(-pi^2*TUBE_RADIUS^2*n[i]^2*pi*TUBE_RADIUS^2*TUBE_LENGTH/PUMP_FLOW/6/Dc)

    # 
    # Ideal cospectra
    f = n[i] * z / U

    if z / L >= 0
        A = 0.124 * (1 + 7.9 * z / L)^0.75
        B = 2.34 * A^-1.1
        UW[i] = f / (A + B * f^2.1)
    else
        if f < 0.24
            UW[i] = 20.78 * f / (1 + 31 * f)^1.575
        else
            UW[i] = 12.66 * f / (1 + 9.6 * f)^2.4
        end
    end

    if z / L >= 0
        A = 0.284 * (1 + 6.4 * z / L)^0.75
        B = 2.34 * A^-1.1
        WX[i] = f / (A + B * f^2.1)
    else
        if f < 0.54
            WX[i] = 12.92 * f / (1 + 26.7 * f)^1.375
        else
            WX[i] = 4.378 * f / (1 + 3.8 * f)^2.4
        end
    end

    # Simpson’s rule weights
    if i == 1
        simp[i] = 1
    elseif iseven(i)
        simp[i] = 4
    else
        simp[i] = 2
    end

    # Final transfer function values
    Suw[i] = simp[i] * sqrt(Twu[i] * Tww[i]) * UW[i]
    Swt[i] = simp[i] * sqrt(Tww[i] * Twt[i]) * WX[i]
    Swq[i] = simp[i] * Ts[i] * sqrt(Tww[i] * Twq[i] * Tdq[i] * Ttq[i]) * WX[i]
    Swc[i] = simp[i] * Ts[i] * sqrt(Tww[i] * Twc[i] * Tdc[i] * Ttc[i]) * WX[i]

    # Prepare next frequency
    if i < 18
        n[i+1] = n[i] * LF
    end
end

# Simpson’s rule integration
TSuw = sum(Suw)
TSwt = sum(Swt)
TSwq = sum(Swq)
TSwc = sum(Swc)

# Correction factors
factor_uw = C / TSuw
factor_wTs = C / TSwt
factor_wq = C / TSwq
factor_wc = C / TSwc

return factor_uw, factor_wTs, factor_wq, factor_wc

end