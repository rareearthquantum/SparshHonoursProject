using DrWatson
@quickactivate "SparshHonoursProject"

include(srcdir("prototype_fft_schemes.jl"))
include(srcdir("plotting.jl"))
include(srcdir("simulation_config.jl"))

default_config = scriptsdir("input_params_large_single_pulse.toml")
config_path = isempty(ARGS) ? default_config : abspath(ARGS[1])
config = load_simulation_config(config_path)

Nz, Nd, Nt, Ny = config.N
Z_range, T_range, Y_range = config.ranges
pulse_params = config.pulse_params
alpha, beta = config.alpha, config.beta

@show config_path config.mode config.name alpha beta

# Run the dimensionless numerical problem.
@time (new_P_grid, new_a_grid),
(new_z_grid, new_t_grid, new_y_grid) =
    evolve_diff_2d_test(
        pulse_params,
        config.N,
        config.ranges,
        (alpha, beta)
    )

plot_z_grid, plot_t_grid, plot_y_grid =
    output_grids(config, (new_z_grid, new_t_grid, new_y_grid))

# Transform the field back to transverse position space.
new_v_grid = fftshift(fftfreq(Ny, 1 / step(new_y_grid)))
for i in eachindex(new_z_grid)
    for j in eachindex(new_t_grid)
        new_a_grid[i, j, :] .*= factor.(new_z_grid[i], new_v_grid, -beta)
    end
end
ifft!(new_a_grid, 3)

# Generate plots. The field remains dimensionless; physical configs restore SI
# coordinates on the axes.
config_tag = replace(config.name, r"[^A-Za-z0-9_-]" => "_")
plot_prefix = plotsdir(config_tag * "_Nz=$(Nz)_Nd=$(Nd)_Nt=$(Nt)_Ny=$(Ny)")

h_numeric = heatmap_myway(
    plot_y_grid,
    plot_z_grid,
    abs2.(new_a_grid[:, end÷2, :]),
    xlabel=config.axis_labels.y,
    ylabel=config.axis_labels.z,
    hlabel="Dimensionless electric-field envelope intensity",
    size=(1000, 1000)
)
savefig(plot_prefix * "_numeric.png")

function analytic_function(z)
    result = zeros(ComplexF64, Nt, Ny)
    for pulse_params_i in pulse_params
        time_params, y_params = pulse_params_i
        propagated_width = sqrt(y_params.width^2 + im * 0.2 * pi^2 * beta * z)
        result .+= pulse.(new_t_grid, time_params...) .* pulse.(
            new_y_grid,
            y_params.center,
            propagated_width,
            y_params.area
        )'
    end
    return result
end

analytic_grid = Array{ComplexF64}(undef, Nz, Nt, Ny)
for i in eachindex(new_z_grid)
    analytic_grid[i, :, :] .= analytic_function(new_z_grid[i])
end

lineplot_inputface = plot(
    plot_y_grid,
    abs2.(analytic_grid[begin, end÷2, :]),
    c=:red,
    label="analytic",
    title="Envelope intensity at input face"
)
plot!(plot_y_grid, abs2.(new_a_grid[begin, end÷2, :]), c=:blue, label="numeric")

lineplot_middle = plot(
    plot_y_grid,
    abs2.(analytic_grid[end÷2, end÷2, :]),
    c=:red,
    label="analytic",
    title="Envelope intensity in the middle of medium"
)
plot!(plot_y_grid, abs2.(new_a_grid[end÷2, end÷2, :]), c=:blue, label="numeric")

lineplot_exitface = plot(
    plot_y_grid,
    abs2.(analytic_grid[end, end÷2, :]),
    c=:red,
    label="analytic",
    title="Envelope intensity at exit face"
)
plot!(plot_y_grid, abs2.(new_a_grid[end, end÷2, :]), c=:blue, label="numeric")

h_analytic = heatmap_myway(
    plot_y_grid,
    plot_z_grid,
    abs2.(analytic_grid[:, end÷2, :]),
    xlabel=config.axis_labels.y,
    ylabel=config.axis_labels.z,
    hlabel="Diffraction-only analytic envelope intensity"
)

plot(h_analytic, h_numeric, layout=(2, 1))
savefig(plot_prefix * "_analytic_comparison.png")

#=
plot(lineplot_inputface, lineplot_middle, lineplot_exitface, layout=(3, 1))
savefig(plot_prefix * "_analytic_lines.png")
=#
