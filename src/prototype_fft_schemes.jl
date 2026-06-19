using OrdinaryDiffEq
using FFTW
using Interpolations

include(srcdir("stepping_schemes.jl"))
include(srcdir("initialise_params.jl"))
include(srcdir("equations.jl"))
include(srcdir("input_pulse_methods.jl"))


function initialise_params_test(
    pulse_params,
    N::NTuple{4,Int64},
    B::NTuple{3,NTuple{2,Float64}}
)

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

    P_grid = Array{ComplexF64}(undef, Nz, Nt, Ny)
    P_grid[:, 1, :] .= zeros(Nz, Ny)

    a_grid = Array{ComplexF64}(undef, Nz, Nt, Ny)
    a_grid[1, :, :] .= initialise_pulses(t_grid, y_grid, pulse_params)
    a_grid[1,:,:] .= fftshift(fft(a_grid[1,:,:],2),2)
    @show size(a_grid)

    grids = (z_grid, t_grid, d_grid, y_grid, P_grid, a_grid)

    return grids
end


function f_a_fft_test!(alpha, P, factor)
    im * alpha * P * factor
end

function factor(z, v, beta)
    return exp(im * beta * 4 * pi^2 * v^2 * z)
end

function factor_grid(z_grid, v_grid, beta)
    return exp.(im .* beta .* 4 .* pi^2 .* v_grid' .^ 2 .* z_grid)
end

function evolve_a_fft_grid_test(t_grid::AbstractVector,
    a_grid::AbstractArray,
    l::Int,
    alpha::Real,
    P_grid::AbstractArray,
    dz::Real,
    i::Int,
    phasefactor_grid::AbstractArray
)

    for k in eachindex(t_grid)
        f(iz) = f_a_fft_test!(alpha, P_grid[iz, k, l], phasefactor_grid[iz, l])
        stepping_ab_2step!(view(a_grid, :, k, l), f, dz, i)
    end

end

function evolve_diff_2d_test(
    pulse_params::AbstractArray,
    N::NTuple{4,Int64},
    B::NTuple{3,NTuple{2,Float64}},
    p::NTuple{2,Float64}
)

    (Ti, Tf) = B[2]
    alpha::Real, beta::Real = p

    z_grid, t_grid, d_grid, y_grid, P_grid, a_grid = initialise_params_test(pulse_params, N, B)

    v_grid::LinRange{Float64,Int64} = fftshift(fftfreq(length(y_grid), 1 / step(y_grid)))
    phasefactor_grid = factor_grid(z_grid, v_grid, beta)
    s0 = zeros(ComplexF64, Nd)

    a_interp_grid = Array{AbstractInterpolation}(undef, Nz, Ny)
    for i in eachindex(z_grid)
        for l in eachindex(y_grid)
            a_interp_grid[i, l] = LinearInterpolation(t_grid, a_grid[i, :, l])
        end
    end

    initial_p = (d_grid[1, :, 1], a_interp_grid[1, 1])
    prob = ODEProblem(s_evolve!, s0, (Ti, Tf), initial_p)
    integrator = init(prob, Tsit5(), reltol=1e-3, abstol=1e-6, saveat=t_grid)

    for i in eachindex(z_grid)
        for l in eachindex(y_grid)

            new_p = (d_grid[i, :, l], a_interp_grid[i, l])
            integrator.p = new_p
            reinit!(integrator, s0)
            solve!(integrator)

            #evolving polarisation
            P_grid[i, :, l] .= dropdims(sum(integrator.sol, dims=1), dims=1)

            #evolving electric field
            evolve_a_fft_grid_test(t_grid, a_grid, l, alpha, P_grid, step(z_grid), i, phasefactor_grid)

        end
    end

    return (P_grid, a_grid), (z_grid, t_grid, y_grid)
end