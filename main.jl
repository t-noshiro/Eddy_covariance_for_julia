###########################################################
#This program was made for EC calculation in julia.
#data will be calculated at your preferred time period.
#Author: Docopypl (git)
###########################################################
#using Revise
#using Profile
@time begin
println("Starting ECsystem...")
include("./src/ECsystem.jl")

import .ECsystem
#print(Threads.nthreads())
ECsystem.main()
#Profile.print()
end