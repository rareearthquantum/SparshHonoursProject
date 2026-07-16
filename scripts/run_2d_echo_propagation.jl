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
parameter_info = "_Nt=$(cfg.Nt)_Nd=$(cfg.Nd)_Nz=$(cfg.Nz)_Ny=$(cfg.Ny)_dwidth=$(cfg.d_width)_"

fig = plot_2d_echo_snapshots(result)
plot_output_dir = joinpath(dirname(@__DIR__), "plots", "prototype_echo_2d")
mkpath(plot_output_dir)
plot_path = joinpath(plot_output_dir, "$(timestamp)" * parameter_info * ".png")
savefig(fig, plot_path)
println("Saved $plot_path")

anim_output_dir = joinpath(plot_output_dir, "anim")
mkpath(anim_output_dir)
anim_path = joinpath(anim_output_dir, "$(timestamp)" * parameter_info * ".gif")
animate_2d_echo(result; filename=anim_path);
println("Saved $anim_path")

data_output_dir = joinpath(dirname(@__DIR__), "data", "prototype_echo_2d")
mkpath(data_output_dir)
data_path = joinpath(data_output_dir, "$(timestamp)" * parameter_info * ".jld2")
save_result(data_path, result, elapsed)
println("Saved $data_path")


Omega_abs2_zmaxes = zeros(Float64,length(result.z_vec))
for i in eachindex(result.z_vec)
    Omega_abs2_zmaxes[i] = maximum(abs2.(result.Omega[:,i,end÷2]))
end

fig2 = plot(result.z_vec, Omega_abs2_zmaxes)
otherplot_output_dir = joinpath(plot_output_dir, "otherplots")
mkpath(otherplot_output_dir)
otherplot_path = joinpath(otherplot_output_dir, "$(timestamp)" * parameter_info * ".png")
savefig(fig2, otherplot_path)
println("Saved $otherplot_path")

display(fig)