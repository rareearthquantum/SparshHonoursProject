using OrdinaryDiffEq
using FFTW
using Interpolations

include(srcdir("stepping_schemes.jl"))
include(srcdir("initialise_params.jl"))
include(srcdir("equations.jl"))


function initialise_params_test(
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
    d_width = 1 / step(t_grid)
    for i in 1:Nz
        for l in 1:Ny
            d_grid[i, :, l] .= LinRange(-d_width, d_width, Nd)
        end
    end

    s_grid = Array{ComplexF64}(undef, Nz, Nd, Nt, Ny)
    s_grid[:, :, 1, :] .= zeros(Nz, Nd, Ny)

    a_grid = Array{ComplexF64}(undef, Nz, Nt, Ny)
    a_grid[1, :, :] .= f_E_t.(t_grid) * fftshift(fft(f_E_y.(y_grid)))'
    @show size(a_grid)

    grids = (z_grid, t_grid, d_grid, y_grid, s_grid, a_grid)

    return grids
end


function f_a_fft_test!(alpha, P, factor)
    im * alpha * P * factor
end

function factor(z, v, beta)
    return exp(im * beta * 4 * pi^2 * v^2 * z)
end

function evolve_a_fft_grid_test(z_grid::LinRange{Float64,Int64},
    v_grid::LinRange{Float64,Int64},
    t_grid::LinRange{Float64,Int64},
    a_grid::Array{ComplexF64},
    l::Int64,
    alpha::Float64,
    s_grid::Array{ComplexF64},
    beta::Float64,
    dz::Float64,
    i::Int64
)

    P = dropdims(sum(s_grid, dims=2), dims=2)

    for k in eachindex(t_grid)
        f(iz) = f_a_fft_test!(alpha, P[iz, k, l], factor(z_grid[iz], v_grid[l], beta))
        a_temp = @view a_grid[:, k, l]
        stepping_ab_2step!(a_temp, f, dz, i)
    end
    #a_grid[i, :, l] .= a_grid[i, :, l] .* factor(z_grid[i], v_grid[l], -beta)
    #a_grid[i, :, l] .= fftshift(ifft(a_grid[i, :, l]))

end

function evolve_diff_2d_test(
    f::NTuple{2,Function},
    N::NTuple{4,Int64},
    B::NTuple{3,NTuple{2,Float64}},
    p::NTuple{2,Float64}
)

    (Ti, Tf) = B[2]

    z_grid, t_grid, d_grid, y_grid, s_grid, a_grid = initialise_params_test(f, N, B)
    v_grid::LinRange{Float64,Int64} = fftshift(fftfreq(Ny, 1 / step(y_grid)))

    alpha::Float64, beta::Float64 = p

    for i in eachindex(z_grid)
        for l in eachindex(y_grid)

            #evolving polarisation
            prob = ODEProblem(s_evolve!, s_grid[i, :, 1, l], (Ti, Tf), (d_grid[i, :, l], LinearInterpolation(t_grid, a_grid[i, :, l])))

            s_grid[i, :, :, l] .= solve(prob, Tsit5(), reltol=1e-3, abstol=1e-6, saveat=t_grid)

            #evolving electric field
            evolve_a_fft_grid_test(z_grid, v_grid, t_grid, a_grid, l, alpha, s_grid, beta, step(z_grid), i)

        end
    end

    return (s_grid, a_grid), (z_grid, t_grid, y_grid)
end