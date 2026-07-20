import Pkg
Pkg.activate(dirname(@__DIR__))

using Dates

include("../src/echo_2d_propagation.jl")
include("../src/echo_2d_plotting.jl")
include("../src/echo_saving.jl")

cfg = EchoConfig()

stats = @timed run_2d_propagation(cfg)
result = stats.value
elapsed = stats.time
println("Took $elapsed seconds to run.")

#saving
timestamp = Dates.format(now(), dateformat"yyyymmdd-HHMMSS-sss")
parameter_info = "_Nt=$(cfg.Nt)_Nd=$(cfg.Nd)_Nz=$(cfg.Nz)_Ny=$(cfg.Ny)_dwidth=$(cfg.d_width)_alpha=$(cfg.alpha)_beta=$(cfg.beta)_pulsecount=$(length(cfg.pulses))_"
for pulse in cfg.pulses
    global parameter_info *= "_area=$(pulse[1].area)_width=$(pulse[1].width)_"
end

fig = plot_2d_echo_snapshots(result)
plot_output_dir = joinpath(dirname(@__DIR__), "plots", "prototype_echo_2d")
mkpath(plot_output_dir)
plot_path = joinpath(plot_output_dir, "$(timestamp)" * parameter_info * ".png")
savefig(fig, plot_path)
println("Saved plot to /plots/")

anim_output_dir = joinpath(plot_output_dir, "anim")
mkpath(anim_output_dir)
anim_path = joinpath(anim_output_dir, "$(timestamp)" * parameter_info * ".gif")
animate_2d_echo(result; filename=anim_path);
println("Saved animation to /plots/anim/")

data_output_dir = joinpath(dirname(@__DIR__), "data", "prototype_echo_2d")
mkpath(data_output_dir)
data_path = joinpath(data_output_dir, "$(timestamp)" * parameter_info * ".jld2")
save_result(data_path, result, elapsed)
println("Saved jld2 data to /data/prototype_echo_2d/")


Omega_abs2_zmaxes = zeros(Float64,length(result.z_vec))
for i in eachindex(result.z_vec)
    Omega_abs2_zmaxes[i] = maximum(abs2.(result.Omega[:,i,end÷2]))
end

fig2 = plot(result.z_vec, Omega_abs2_zmaxes, xlabel="z", ylabel="time max of intensity at y=0", ylims=(0.8*minimum(Omega_abs2_zmaxes),1.2*maximum(Omega_abs2_zmaxes)))
otherplot_output_dir = joinpath(plot_output_dir, "otherplots")
mkpath(otherplot_output_dir)
otherplot_path = joinpath(otherplot_output_dir, "time_max_vs_z_$(timestamp)" * parameter_info * ".png")
savefig(fig2, otherplot_path)
println("Saved another plot to /plots/otherplots/")

fig_energy_transmission = plot_total_energy_vs_optical_depth(result)
energy_transmission_plot_path = joinpath(otherplot_output_dir, "energy_transmission_$(timestamp)" * parameter_info * ".png")
savefig(fig_energy_transmission, energy_transmission_plot_path)
println("Saved energy transmission plot to /plots/otherplots/")