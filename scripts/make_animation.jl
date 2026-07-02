using DrWatson
@quickactivate "SparshHonoursProject"

include(srcdir("input_pulse_methods.jl"))
include(srcdir("echo_propagation.jl"))

cfg = EchoConfig()

@time result = run_propagation(cfg)

animate_field_and_polarisation(result; filename="field_and_polarisation.gif", fps=30)
