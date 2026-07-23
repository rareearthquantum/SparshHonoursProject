import Pkg
Pkg.activate(dirname(@__DIR__))

using Dates

include("../src/echo_2d_propagation.jl")
include("../src/echo_2d_plotting.jl")
include("../src/echo_saving.jl")

# Init config
cfg = EchoConfig()

run_simulation(cfg)