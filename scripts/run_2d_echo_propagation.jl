import Pkg
Pkg.activate(dirname(@__DIR__))

using Dates

include("../src/echo_2d_propagation.jl")
include("../src/echo_2d_plotting.jl")
include("../src/echo_saving.jl")

# Init config
cfg = EchoConfig()

# RUN
stats = @timed run_2d_propagation(cfg)

# Take results
result = stats.value
elapsed = stats.time
println("Took $elapsed seconds to run.")

# Setting info for saving into file
timestamp = Dates.format(now(), dateformat"yyyymmdd-HHMMSS-sss")
parameter_info = "_Nt=$(cfg.Nt)_Nd=$(cfg.Nd)_Nz=$(cfg.Nz)_Ny=$(cfg.Ny)_dwidth=$(cfg.d_width)_alpha=$(cfg.alpha)_beta=$(cfg.beta)_pulsecount=$(length(cfg.pulses))_"
for pulse in cfg.pulses
    global parameter_info *= "_area=$(pulse[1].area)_width=$(pulse[1].width)_"
end


# Setting plot directory and helper function
plot_output_dir = joinpath(dirname(@__DIR__), "plots", "prop_2d")
mkpath(plot_output_dir)
plot_n_save(func, name) = save_plot(result, func, plot_output_dir, name; parameter_info=parameter_info, timestamp=timestamp)


# Plotting
#plot_n_save(plot_sum_omega, "energy")
#plot_n_save(x -> plot_sum_omega(x; operation=real), "area")
plot_n_save(x -> plot_soliton_lineshapes(x;nslices=5), "soliton_lineshapes")
#plot_n_save(animate_field_2d, "anim")

# Saving jld2 data
#save_data(result, elapsed, "prop_2d")