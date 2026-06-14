using DrWatson
@quickactivate "SparshHonoursProject"

include(srcdir("prototype_fft_schemes.jl"))
include(srcdir("input_pulse_methods.jl"))
include(srcdir("plotting.jl"))

##INPUT PARAMS

Nz = 128 #number of z-sites
Nd = 2 #number of detunings/atoms at each z-site
Nt = 32 #number of time steps
Ny = 32

Z_range = (0.0, 1.0)
T_range = (0.0, 10.0)
Y_range = (-5.0, 5.0)

alpha = 0.0
beta = 3.0

@show "test run"

##Run
@time (new_s_grid, new_a_grid),
(new_z_grid, new_t_grid, new_y_grid) =
    evolve_diff_2d_test(
        (Ein_t, Ein_y),
        (Nz, Nd, Nt, Ny),
        (Z_range, T_range, Y_range),
        (alpha, beta)
    )


##IFFT back
new_v_grid = fftshift(fftfreq(Ny, 1 / step(new_y_grid)))


for i in eachindex(new_z_grid)
    for j in eachindex(new_t_grid)
        new_a_grid[i, j, :] .= new_a_grid[i, j, :] .* factor.(new_z_grid[i], new_v_grid, -beta)
        new_a_grid[i, j, :] .= ifft(new_a_grid[i, j, :])
    end
end




##Generate plots
plotting_a_intensity(new_t_grid, new_y_grid, new_z_grid, new_a_grid,clims=(0.0,1.0))
savefig("plots/test_a_abs2" *
        "_alpha=" * string(alpha) *
        "_beta=" * string(beta) *
        "_Nz=" * string(Nz) *
        "_Nd=" * string(Nd) *
        "_Nt=" * string(Nt) *
        "_Ny=" * string(Ny) *
        ".png")