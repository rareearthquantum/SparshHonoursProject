import Pkg
Pkg.activate(dirname(@__DIR__))

using Dates

include("../src/echo_2d_propagation.jl")
include("../src/echo_2d_plotting.jl")

cfg = EchoConfig()

@time result = run_2d_propagation(cfg)

fig = plot_2d_echo_snapshots(result)

timestamp = Dates.format(now(), dateformat"yyyymmdd-HHMMSS-sss")
output_dir = joinpath(dirname(@__DIR__), "plots", "prototype_echo_2d")
mkpath(output_dir)
plot_path = joinpath(output_dir, "$(timestamp)" * ".png")
savefig(fig, plot_path)

display(fig)
sleep(5)

plot_path = joinpath(output_dir, "$(timestamp)" * ".gif")
animate_2d_echo(result; filename=plot_path)