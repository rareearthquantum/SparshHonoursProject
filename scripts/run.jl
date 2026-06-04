using DrWatson
@quickactivate "SparshHonoursProject"

##INPUT PARAMS

Nz = 128 * 2 #number of z-sites
Nd = 2 #number of detunings/atoms at each z-site
Nt = 16 #number of time steps
Ny = 32

Zi, Zf = 0.0, 10.0
Ti, Tf = 2.5, 7.5
Yi, Yf = -5.0, 5.0

alpha = 0.0
beta = 5.0e-2

##Run

(new_s_grid, new_a_grid),
(new_z_grid, new_t_grid, new_y_grid) =
    evolve_diff_2d(
        (Ein_t, Ein_y),
        (Nz, Nd, Nt, Ny),
        (
            (Zi, Zf),
            (Ti, Tf),
            (Yi, Yf)
        ),
        (alpha, beta)
    );

