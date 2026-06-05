function initialise_params(
    f::NTuple{2,Function},
    N::NTuple{4,Int64},
    B::NTuple{3,NTuple{2,Float64}}
)

    f_E_t, f_E_y = f
    Nz, Nd, Nt, Ny = N
    (Zi, Zf), (Ti, Tf), (Yi, Yf) = B

    ### Grids
    ## Z
    z_grid = LinRange(Zi, Zf, Nz)
    dz = (Zf - Zi) / Nz

    ## Y
    y_grid = LinRange(Yi, Yf, Ny)
    dy = (Yf - Yi) / Ny

    ## T
    t_grid = LinRange(Ti, Tf, Nt)
    dt = (Tf - Ti) / Nt

    ## D
    d_grid = Array{Float64}(undef, Nz, Nd, Ny)
    df = 1 / (Nt * dt)
    d_width = df * Nt
    for i in 1:Nz
        for l in 1:Ny
            d_grid[i, :, l] = LinRange(-d_width, d_width, Nd)
        end
    end


    ## S
    s_grid = Array{ComplexF64}(undef, Nz, Nd, Nt, Ny)
    for l in 1:Ny
        s_grid[:, :, 1, l] .= zeros(Nz, Nd)
    end

    ## A
    a_grid = Array{ComplexF64}(undef, Nz, Nt, Ny)
    for l in 1:Ny
        a_grid[1, :, l] .= f_E_t.(t_grid) .* f_E_y(y_grid[l])
    end

    ## V
    v_freq = fftfreq(Ny, 1 / dy)
    v_factor = -4.0 .* π^2 .* v_freq .^ 2


    grids = (z_grid, t_grid, d_grid, y_grid, s_grid, a_grid)
    differentials = (dz, dt, dy)
    v = (v_freq, v_factor)

    return (grids, differentials, v)
end
