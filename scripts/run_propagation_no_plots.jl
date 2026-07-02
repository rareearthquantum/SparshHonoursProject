using DrWatson
@quickactivate "SparshHonoursProject"

include(srcdir("input_pulse_methods.jl"))
include(srcdir("echo_propagation.jl"))

cfg = EchoConfig()

@time result = run_propagation(cfg)

# Useful variables after include("scripts/run_propagation_no_plots.jl"):
# result.time_vec, result.z_vec, result.Omega, result.P, result.detunings
