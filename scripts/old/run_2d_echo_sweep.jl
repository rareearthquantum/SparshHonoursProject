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
parameter_info = "_Nt=$(cfg.Nt)_Nd=$(cfg.Nd)_Nz=$(cfg.Nz)_Ny=$(cfg.Ny)_dwidth=$(cfg.d_width)_"
plot_path = joinpath(output_dir, "$(timestamp)" * parameter_info * ".png")
savefig(fig, plot_path)

plot_path = joinpath(output_dir, "anim", "$(timestamp)" * parameter_info * ".gif")
animate_2d_echo(result; filename=plot_path)

display(fig)