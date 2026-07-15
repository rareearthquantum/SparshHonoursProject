import Pkg
Pkg.activate(dirname(@__DIR__))

using Dates

include("../src/echo_propagation.jl")

cfg = EchoConfig()

@time result = run_propagation(cfg)

fig = plot_propagation_super_compact(result)

timestamp = Dates.format(now(), dateformat"yyyymmdd-HHMMSS-sss")
output_dir = joinpath(dirname(@__DIR__), "plots", "echo_1d_propagation")
mkpath(output_dir)
plot_path = joinpath(output_dir, "$(timestamp)" * ".png")
savefig(fig, plot_path)
println("Saved $plot_path")

if (length(cfg.pulses) > 1)
    fig_echo_transmission_and_efficiency = plot_echo_transmission_and_efficiency_vs_optical_depth(result)

    timestamp = Dates.format(now(), dateformat"yyyymmdd-HHMMSS-sss")
    plot_path = joinpath(output_dir, "$(timestamp)" * ".png")
    savefig(fig_echo_transmission_and_efficiency, plot_path)
    println("Saved $plot_path")
end


display(fig)