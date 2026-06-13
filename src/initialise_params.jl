using FFTW

function initialise_params(
    f::NTuple{2,Function},
    N::NTuple{4,Int64},
    B::NTuple{3,NTuple{2,Float64}}
)

    f_E_t, f_E_y = f
    Nz, Nd, Nt, Ny = N
    (Zi, Zf), (Ti, Tf), (Yi, Yf) = B



    z_grid = LinRange(Zi, Zf, Nz)
    y_grid = LinRange(Yi, Yf, Ny)
    t_grid = LinRange(Ti, Tf, Nt)

    d_grid = Array{Float64}(undef, Nz, Nd, Ny)
    d_width = 1 / (step(t_grid))
    for i in 1:Nz
        for l in 1:Ny
            d_grid[i, :, l] .= LinRange(-d_width, d_width, Nd)
        end
    end

    s_grid = Array{ComplexF64}(undef, Nz, Nd, Nt, Ny)
    #setting boundary conditions
    s_grid[:, :, 1, :] .= zeros(Nz, Nd, Ny)

    a_grid = Array{ComplexF64}(undef, Nz, Nt, Ny)
    #setting boundary conditions
    a_grid[1, :, :] .= f_E_t.(t_grid) * (f_E_y.(y_grid))'

    v_freq = fftfreq(Ny, 1 / step(y_grid))
    v_factor = -4.0 .* π^2 .* v_freq .^ 2

    grids = (z_grid, t_grid, d_grid, y_grid, s_grid, a_grid)
    v_stuff = (v_freq, v_factor)

    return grids, v_stuff
end