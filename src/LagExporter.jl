#module LagExporter

#using ..Config

#export LAG_MAX, LAG_MIN, LAG_OFFSET, LAG_OFFSET_ANOTHER

#include(joinpath(@__DIR__,"config_module.jl"))

const LAG_MAX = LAG_MAX_SEC * SAMPLING_SPAN
const LAG_MIN = LAG_MIN_SEC * SAMPLING_SPAN
const LAG_OFFSET = LAG_OFFSET_SEC * SAMPLING_SPAN
const LAG_OFFSET_ANOTHER = LAG_OFFSET_ANOTHER_SEC * SAMPLING_SPAN

#end