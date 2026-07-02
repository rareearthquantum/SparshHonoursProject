# Main include file for the propagation code.
#
# In a DrWatson project, load this with:
#     include(srcdir("echo_propagation.jl"))
# after loading input_pulse_methods.jl.

include("echo_config.jl")
include("echo_equations.jl")
include("echo_integrators.jl")
include("echo_simulation.jl")
include("echo_plotting.jl")
