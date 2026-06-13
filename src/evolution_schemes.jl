using OrdinaryDiffEq
using FFTW
using Interpolations

include(srcdir("stepping_schemes.jl"))
include(srcdir("initialise_params.jl"))
include(srcdir("equations.jl"))

function evolve_a_grid(
    t_grid::LinRange{Float64,Int64},
    a_grid::Array{ComplexF64},
    v_factor,
    l::Int64,
    alpha::Float64,
    s_grid::Array{ComplexF64},
    beta::Float64,
    dz::Float64,
    i::Int64,
    Nz::Int64
)
    @inbounds for k in eachindex(t_grid)
        dy2 = x -> (ifft(v_factor .* fft(a_grid[x, k, :])))[l]
        f_2d = x -> f_a_2d!(alpha, s_grid[x, :, k, l], beta, dy2(x))
        a_temp = @view a_grid[:, k, l]
        stepping_ab_2step!(a_temp, f_2d, dz, i)
    end
end


function evolve_diff_2d(
    f::NTuple{2,Function},
    N::NTuple{4,Int64},
    B::NTuple{3,NTuple{2,Float64}},
    p::NTuple{2,Float64}
)

    Nz, Nd, Nt, Ny = N
    (Ti, Tf) = B[2]

    (z_grid, t_grid, d_grid, y_grid, s_grid, a_grid), (v_freq, v_factor) = initialise_params(f, N, B)

    alpha::Float64, beta::Float64 = p

    for i in eachindex(z_grid)
        for l in eachindex(y_grid)

            #evolving polarisation
            prob = ODEProblem(s_evolve!, s_grid[i, :, 1, l], (Ti, Tf), (d_grid[i, :, l], LinearInterpolation(t_grid, a_grid[i, :, l])))

            s_grid[i, :, :, l] .= solve(prob, Tsit5(), reltol=1e-3, abstol=1e-6, saveat=t_grid)

            #evolving electric field
            evolve_a_grid(t_grid, a_grid, v_factor, l, alpha, s_grid, beta, step(z_grid), i, Nz)

        end
    end

    return (s_grid, a_grid), (z_grid, t_grid, y_grid)
end