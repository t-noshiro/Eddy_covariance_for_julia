module Config
include(joinpath(@__DIR__,"..","config.jl"))
for name in names(@__MODULE__, all = true)
    if name != :Config && !startswith(string(name), "#")
        @eval export $name
    end 
end
end