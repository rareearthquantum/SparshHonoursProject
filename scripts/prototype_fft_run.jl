using DrWatson
@quickactivate "SparshHonoursProject"

include(srcdir("prototype_fft_schemes.jl"))
include(srcdir("plotting.jl"))
include(srcdir("constants.jl"))


#INPUT PARAMS

const Nz = 64 #number of z-sites
const Nd = 64  #number of detunings/atoms at each z-site
const Nt = 64  #number of time steps
const Ny = 64 *2 

const Z_length = 2milli
const Z_range = (0.0, Z_length)
const T_range = (0.0, 50micro)

const y_pulse_width = 10micro
const ZR = 0.5 * k * y_pulse_width^2
const Y_width = real(sqrt(y_pulse_width^2 + im * 2 * beta * Z_length))
const Y_range = (-Y_width / 2, Y_width / 2)

const seperation = 40micro
const pulse_params = [
    ((center=T_range[2] / 2, width=10micro, area=1.0), (center=0.0, width=y_pulse_width, area=1.0))
]

@show "test run"

#Run
@time (new_P_grid, new_a_grid),
(new_z_grid, new_t_grid, new_y_grid) =
    evolve_diff_2d_test(
        pulse_params,
        (Nz, Nd, Nt, Ny),
        (Z_range, T_range, Y_range),
        (alpha, beta)
    )


# ifft back
new_v_grid = fftshift(fftfreq(Ny, 1 / step(new_y_grid)))
for i in eachindex(new_z_grid)
    for j in eachindex(new_t_grid)
        new_a_grid[i, j, :] .*= factor.(new_z_grid[i], new_v_grid, -beta)
    end
end
ifft!(new_a_grid, 3)


##Generate plots
h_numeric = heatmap_myway(new_y_grid, new_z_grid, abs2.(new_a_grid[:, end÷2, :]),
    xlabel="y", ylabel="z", hlabel="Electric field envelope intensity", size=(1000, 1000))
savefig("plots/test_a_abs2" *
        "_Nz=" * string(Nz) *
        "_Nd=" * string(Nd) *
        "_Nt=" * string(Nt) *
        "_Ny=" * string(Ny) *
        "_Zw=" * string(Z_range[2] - Z_range[1]) *
        "_Tw=" * string(T_range[2] - T_range[1]) *
        "_Yw=" * string(Y_range[2] - Y_range[1]) *
        ".png")

analytic_width = z -> sqrt(y_pulse_width^2 + im * 2 * beta * z)
analytic_pulse_params = z -> ((center=(T_range[2] - T_range[1]) / 2, width=10micro, area=1.0), (center=0.0, width=analytic_width(z), area=1.0))
analytic_function = z -> pulse.(new_t_grid, analytic_pulse_params(z)[1]...) .* pulse.(new_y_grid, analytic_pulse_params(z)[2]...)'

analytic_grid = Array{Complex}(undef, Nz, Nt, Ny)
for i in eachindex(new_z_grid)
    analytic_grid[i, :, :] .= analytic_function(new_z_grid[i])
end

lineplot_inputface = lineplot_exitface = plot(new_y_grid, abs2.(analytic_grid[begin, end÷2, :]), c=:red, label="analytic", title="Envelope intensity at input face")
plot!(new_y_grid, abs2.(new_a_grid[begin, end÷2, :]), c=:blue, label="numeric")

lineplot_middle = plot(new_y_grid, abs2.(analytic_grid[end÷2, end÷2, :]), c=:red, label="analytic", title="Envelope intensity in the middle of medium")
plot!(new_y_grid, abs2.(new_a_grid[end÷2, end÷2, :]), c=:blue, label="numeric")

lineplot_exitface = plot(new_y_grid, abs2.(analytic_grid[end, end÷2, :]), c=:red, label="analytic", title="Envelope intensity at exit face")
plot!(new_y_grid, abs2.(new_a_grid[end, end÷2, :]), c=:blue, label="numeric")

h_analytic = heatmap_myway(new_y_grid, new_z_grid, abs2.(analytic_grid[:, end÷2, :]), xlabel="y", ylabel="z", hlabel="Analytic envelope intensity")

plot(h_analytic, h_numeric, layout=(2, 1))
savefig("plots/analytic__heatmap_abs2_comparision" *
        "_Nz=" * string(Nz) *
        "_Nd=" * string(Nd) *
        "_Nt=" * string(Nt) *
        "_Ny=" * string(Ny) *
        "_Zw=" * string(Z_range[2] - Z_range[1]) *
        "_Tw=" * string(T_range[2] - T_range[1]) *
        "_Yw=" * string(Y_range[2] - Y_range[1]) *
        ".png")

plot(lineplot_inputface, lineplot_middle, lineplot_exitface, layout=(3, 1))
savefig("plots/analytic_lines_abs2_comparision" *
        "_Nz=" * string(Nz) *
        "_Nd=" * string(Nd) *
        "_Nt=" * string(Nt) *
        "_Ny=" * string(Ny) *
        "_Zw=" * string(Z_range[2] - Z_range[1]) *
        "_Tw=" * string(T_range[2] - T_range[1]) *
        "_Yw=" * string(Y_range[2] - Y_range[1]) *
        ".png")