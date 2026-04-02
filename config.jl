###########################################################
#This program was made for EC calculating in julia.
#data will be calculated at your preferred time period.
#Author: Docopypl (git)
###########################################################

###################
###set_directory###
###################
#const USER_NAME = "Null"
const SYSTEM_DIRECTORY = @__DIR__

####################
###set_parameter####
####################
const TWO_IRGA_FLAG = true #If you use two IRGAs, Here will be ture. Otherwise, here is false.
const IS_IRGA1_CP = true
const IS_IRGA2_CP = false #->this mean ANOTHER_*_NAME 's IRGA is CP-IRGA (or not).Ignore this if one irga was used.
const SAMPLING_SPAN = 20 ::Int64#Hz
const TIME_SPAN = 30 ::Int64 #min : Time interval for one period in flux calculation
const FILE_KEYWORD = ".dat" #This is endwith_keyword. Input the word by which only target files can be hit. 
const LAG_MAX_SEC = 1 #seconds
const LAG_MIN_SEC = -2 #seconds
const LAG_OFFSET_SEC = 0 #seconds; this means that timeseriesdata will be slided according to lag_offset.
#lag >0 means that xdata delayed  
const LAG_OFFSET_ANOTHER_SEC = 0
    
const DATA_DIRECTORY =joinpath(SYSTEM_DIRECTORY, "rawdata") #Directory of raw data used when calculating the flux
const RESULT_FILE_NAME = "flux_summary"

const TEST_FLAG = 0 #For testing this system's performance; This calculate the flux using the "testdata" folder's data.

const RESULT_OUT_DIRECTORY = joinpath(SYSTEM_DIRECTORY, "out")

const SOS_COLUMN_is = "Ts" #"sos" = SOS_COLUMN is sos(sonic of sound), "Ts" = SOS_COLUMN is t_sonic. Other values are not permitted.

#CSV_line
const HEADER_LINE = 2
const START_LINE = 5

#CSV_colnames
const TIME_COLUMN_NAME = "TIMESTAMP" 
const CO2_DENSITY_COLUMN_NAME = "CO2"
const WATER_VAPOUR_DENSITY_COLUMN_NAME = "H2O"
const U_COLUMN_NAME = "U" #X-Axis wind speed 
const V_COLUMN_NAME = "V" #Y-Axis wind speed
const W_COLUMN_NAME = "W" #Z-Axis wind speed
const SOS_COLUMN_NAME = "Ts" #Speed of sonic or Sonic temperature
const PRESSURE_COLUMN_NAME = "Press" #Cell pressure
const T_IN_COLUMN_NAME = "Tin" #LI7200's Tin
const T_OUT_COLUMN_NAME = "Tout" #LI7200's Tout

const ANOTHER_CO2_DENSITY_COLUMN_NAME = "CO22"
const ANOTHER_WATER_VAPOUR_DENSITY_COLUMN_NAME = "H2O2"
const ANOTHER_PRESSURE_COLUMN_NAME = "Press"
const ANOTHER_T_IN_COLUMN_NAME = "Ts"
const ANOTHER_T_OUT_COLUMN_NAME = "Ts"

#DateFormat
const IS_SEPARATING_DATE_AND_TIME = false
const DATE_COLUMN_NAME = "Date"
const DATE_FORMAT = "yyyy-mm-dd"
const TIME_FORMAT = "HH:MM:SS.sss"
const SEPARATING_CHARACTER = " "
#If only time colmun is used, please set up the variables below as follows:
#DateFormat: yyyy-mm-ddTHH:MM:SS.sss 
#-> DATE_FORMAT = "yyyy-mm-dd"
#-> TIME_FORMAT = "HH:MM:SS.sss"
#-> SEPARATING_CHARACTER = "T"

#Convert_unit_to_SI_standard_units (m, kg, s, K, mol, Pa etc.)
const CONVERT_UNIT_CO2_DENSITY = 1e-3 #to convert CO2 density unit to mol m^-3
const CONVERT_UNIT_WATER_VAPOUR_DENSITY = 1e-3 #to convert H2O density unit to mol m^-3
const CONVERT_UNIT_PRESSURE = 1e3 #to convert pressure unit to Pa

const CONVERT_UNIT_ANOTHER_CO2_DENSITY = 1e-3 #to convert CO2 density unit to mol m^-3
const CONVERT_UNIT_ANOTHER_WATER_VAPOUR_DENSITY = 1e-3 #to convert H2O density unit to mol m^-3
const CONVERT_UNIT_ANOTHER_PRESSURE = 1e3 #to convert pressure unit to Pa

const PRESSURE_ATMOSPHERE = 1.01*1e5 #Pa

const ATMOSPHERIC_PARAMETER_IRGA = 2 #NEITHER DRYING NOR MODIFIED IRGA No.

####################################################################
####################################################################
#Moncrieff97
#IRGA1
const MONCRIEFF97_FLAG_1 = true
const TUBE_RADIUS_1 = 5.33e-3/2 #m
const TUBE_LENGTH_1 = 0.72 #m
const PUMP_FLOW_1 = 1.5e-4 #m3 s-1
const SAT_PATH_LEN = 0.1 #m
const IRGA_PATH_LEN_1 = 0.13 #
const Z_OBSERVE_1 = 1.63 #Observation height of IRGA1
#IRGA2
const MONCRIEFF97_FLAG_2 = false
const TUBE_RADIUS_2 = 5.33e-3/2 #m
const TUBE_LENGTH_2 = 0.72 #m
const PUMP_FLOW_2 = 1.5e-4 #m3 s-1
const IRGA_PATH_LEN_2 = 0.13 #m
const Z_OBSERVE_2 = 1.63 #Observation height of IRGA2

const ADDITIONAL_FACTOR = 1.05 #Noshiro et al., 2024: 1.05 #read this paper to use this factor. Otherwise this should be 1.

const IRGA_SEPARATIONS_1 = [
    1e-2*sqrt(8^2+21^2+1^2)
    1e-2*sqrt(10^2+21^2+1^2)
    1e-2*sqrt(23^2+11^2+1^2)
]

const IRGA_SEPARATIONS_2 = [
    1e-2*sqrt(8^2+21^2+1^2)
    1e-2*sqrt(10^2+21^2+1^2)
    1e-2*sqrt(23^2+11^2+1^2)
]

#format: yyyy-mm-ddTHH:MM:SS <- FOLLOW THIS FORMAT EXACTLY. NO EXCEPTIONS. 
const PERIOD_STARTS = [
    "2023-05-13T12:00:00"
    "2023-05-16T16:00:00"
    "2023-05-24T17:00:00"
]

const PERIOD_ENDS = [
    "2023-05-16T16:00:00"
    "2023-05-24T17:00:00"
    "2023-05-31T23:59:00"
]


#basic constants
const R = 8.314   # gas constant (m2 kg s-2 K-1 mol-1)
const T0 = 273.15  # 0degC in kelvin
const Md = 28.97   # g/mol of dry air
const Mw = 18.015  # g/mol of h2o
# const Mc= 44.01   # g/mol of co2
const Cpd = 1.00467 # heat capacity of dry air [J/(g K)]
const Cpv = 1.875   # heat capacity of water vapor [J/(g K)]

