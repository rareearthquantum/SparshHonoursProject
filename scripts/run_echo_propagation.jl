import Pkg
Pkg.activate(dirname(@__DIR__))

include("../src/input_pulse_methods.jl")
include("../src/echo_propagation.jl")

cfg = EchoConfig()

@time result = run_propagation(cfg)

fig = plot_propagation_compact(result)
display(fig)
