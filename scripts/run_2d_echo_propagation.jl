println("Script running.")
start_time = time()

import Pkg
Pkg.activate(dirname(@__DIR__))

using Dates

include("../src/echo_2d_propagation.jl")

# Init config
cfg = EchoConfig()
println("Config initialised after $(current_runtime(start_time))s.")

run_2d_simulation(cfg)