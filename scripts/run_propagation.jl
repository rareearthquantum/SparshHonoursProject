using DrWatson
@quickactivate "SparshHonoursProject"

include(srcdir("input_pulse_methods.jl"))
include(srcdir("echo_propagation.jl"))

cfg = EchoConfig()

@time result = run_propagation(cfg)

fig = plot_propagation(result)
display(fig)
