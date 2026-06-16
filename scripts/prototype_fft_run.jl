using DrWatson
@quickactivate "SparshHonoursProject"

include(srcdir("prototype_fft_schemes.jl"))
include(srcdir("plotting.jl"))

const frequency = 200e12
const omega0 = 2.0 * pi * frequency
const c = 3e8
const k = omega0 / c


#INPUT PARAMS

const Nz = 128 * 2^2  #number of z-sites
const Nd = 2 #number of detunings/atoms at each z-site
const Nt = 32 * 2^3 #number of time steps
const Ny = 32 * 2^5

const Z_range = (0.0, 10e-2)
const T_range = (0.0, 1.0e-4)
const Y_width = 10e-3
const Y_range = (-Y_width / 2, Y_width / 2)

const alpha = 0.0
const beta = 1 / (2k)

#center,width,area
const seperation = 40e-6
const pulse_params = [
    (PulseParams(T_range[2] / 2, 10e-6, 1.0), PulseParams(-seperation / 2, 10e-6, 1.0e-6)),
    (PulseParams(T_range[2] / 2, 10e-6, 1.0), PulseParams(seperation / 2, 10e-6, 1.0e-6))
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
@time for i in eachindex(new_z_grid)
    for j in eachindex(new_t_grid)
        new_a_grid[i, j, :] .*= factor.(new_z_grid[i], new_v_grid, -beta)
    end
end
ifft!(new_a_grid, 3)


##Generate plots
@time plotting_a_intensity(new_t_grid, new_y_grid, new_z_grid, new_a_grid, clims=(0.0, 1.0e7))
savefig("plots/test_a_abs2" *
        "_alpha=" * string(alpha) *
        "_beta=" * string(beta) *
        "_Nz=" * string(Nz) *
        "_Nd=" * string(Nd) *
        "_Nt=" * string(Nt) *
        "_Ny=" * string(Ny) *
        ".png")





#=
For soon:
- y pulse width needs to be 10 microns, so 10^-6 m
- In time its 10 microseconds, 10^-6 s
- For that I need to be able to input pulse parameters
- Also to get two 10 micron pulses in y seperated by 40 microns
- 10cm propagation in z
- also make analytic plots of everything so far, so: attenuation in 1 spatial dimension, 2 spatial dimensions, diffraction
- also only store P(y,z,t) and evolve single sigma minus for fixed y=y' and z=z', adding the resulting 1 vector to P(y',z',t) for all detunings at y',z'

=#