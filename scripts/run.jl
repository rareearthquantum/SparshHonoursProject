using DrWatson
@quickactivate "SparshHonoursProject"

include(srcdir("evolution_schemes.jl"))
include(srcdir("input_pulse_methods.jl"))
include(srcdir("plotting.jl"))

##INPUT PARAMS

Nz = 64 #number of z-sites
Nd = 2 #number of detunings/atoms at each z-site
Nt = 16 #number of time steps
Ny = 16

Z_range = (0.0, 10.0)
T_range = (0.0, 10.0)
Y_range = (-5.0, 5.0)

alpha = 0.0
beta = 1.0e-2

@show "input params initialised"

##Run
@time (new_s_grid, new_a_grid),
(new_z_grid, new_t_grid, new_y_grid) =
    evolve_diff_2d(
        (Ein_t, Ein_y),
        (Nz, Nd, Nt, Ny),
        (Z_range, T_range, Y_range),
        (alpha, beta)
    )

##Generate plot
plot_a_intensity_vs_y_z(new_y_grid, new_z_grid, new_a_grid)
savefig("plots/myplot.png")