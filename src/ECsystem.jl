#Author Taichi Noshiro
module ECsystem

using Dates, Statistics, DataFrames, CSV, Interpolations

###################
###set_directory###
###################

include("../config.jl")
include("LagExporter.jl")
#basic functions
include("basic_math.jl")
include("dsp_ver4.jl")
#structure and bit complex functions
include("utilities.jl")
include("xcov_ver2.jl")
include("axisrotation_fix.jl")
include("Obukhov_length_ver2.jl")
include("Moncrieff97.jl")
include("winddirection_ver2.jl")

#QC
include("stationarity_test.jl")
include("rannik2016_ver2.jl")

#processor
include("Fluxcal_ver4.jl")
include("makingdata_ver5.0.jl")
include("make_results.jl")
include("output_results.jl")

function main()
    data_dir = DATA_DIRECTORY
    if TEST_FLAG == 1
        data_dir = joinpath(SYSTEM_DIRECTORY, "testdata")
    end
    data_files = search_file(data_dir, FILE_KEYWORD)
    buffer = FluxBuffer_empty()
    make_results(data_files, buffer)
    output_results(buffer)

    
    println("compl.")
end

export main

end # module 