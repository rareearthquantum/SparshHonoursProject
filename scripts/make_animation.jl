import Pkg
Pkg.activate(dirname(@__DIR__))

include("../src/input_pulse_methods.jl")
include("../src/echo_propagation.jl")

cfg = EchoConfig()

@time result = run_propagation(cfg)

animate_field_and_polarisation(result; filename="field_and_polarisation.gif", fps=30)
