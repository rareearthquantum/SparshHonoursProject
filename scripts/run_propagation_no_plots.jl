import Pkg
Pkg.activate(dirname(@__DIR__))

include("../src/echo_propagation.jl")

cfg = EchoConfig()

@time result = run_propagation(cfg)

# Useful variables after include("scripts/run_propagation_no_plots.jl"):
# result.time_vec, result.z_vec, result.Omega, result.P, result.detunings
