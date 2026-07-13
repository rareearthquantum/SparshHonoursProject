import Pkg
Pkg.activate(dirname(@__DIR__))

using Dates

include("../src/echo_propagation.jl")

const DEFAULT_PULSE_PROFILE = joinpath(
    @__DIR__, "parameters", "echo_pulses_auto.jl")

profile_path = isempty(ARGS) ? DEFAULT_PULSE_PROFILE : abspath(first(ARGS))
cfg = include(profile_path)

@show cfg
@time result = run_propagation(cfg)

fig = plot_propagation_compact(result)
display(fig)

#plot
timestamp = Dates.format(now(), dateformat"yyyymmdd-HHMMSS-sss")
output_dir = joinpath(dirname(@__DIR__), "plots", "auto_pulses")
mkpath(output_dir)
plot_path = joinpath(output_dir, "$(timestamp)" * ".png")
savefig(fig, plot_path)
println("Saved $plot_path")