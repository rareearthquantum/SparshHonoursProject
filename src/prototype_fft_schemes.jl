using FFTW

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
    P_grid[:, 1, :] .= zeros(Complex, Nz, Ny)

    a_grid = Array{ComplexF64}(undef, Nz, Nt, Ny)
    a_grid[1, :, :] .= initialise_pulses(t_grid, y_grid, pulse_params)
    a_grid[1, :, :] .= fftshift(fft(a_grid[1, :, :], 2), 2)
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

function evolve_polarisation_rk4!(
    P::AbstractVector,
    d::AbstractVector,
    a::AbstractVector,
    dt::Real
)
    length(P) == length(a) || throw(DimensionMismatch("P and a must have the same length"))
    length(P) > 0 || return P


    s = zeros(ComplexF64, length(d))
    k1 = similar(s)
    k2 = similar(s)
    k3 = similar(s)
    k4 = similar(s)

    d_width = d[end]-d[begin]
    d_step = d_width / length(d)
    weights = ones(Real, Nd) ./ d_width

    for j in 1:(length(a)-1)
        # RK4 evaluates the forcing at the half step.  The adjacent field
        # samples give its second-order midpoint value without an interpolator.
        a_mid = (a[j] + a[j+1]) / 2

        @. k1 = im * (d * s + a[j])
        @. k2 = im * (d * (s + (dt / 2) * k1) + a_mid)
        @. k3 = im * (d * (s + (dt / 2) * k2) + a_mid)
        @. k4 = im * (d * (s + dt * k3) + a[j+1])
        @. s += (dt / 6) * (k1 + 2k2 + 2k3 + k4)

        # Uniform normalized detuning distribution. This keeps the physical
        # coupling independent of the numerical detuning-grid size Nd.
        P[j+1] = sum(weights .* s) * d_step
    end

    return P
end



function evolve_diff_2d_test(
    pulse_params::AbstractArray,
    N::NTuple{4,Int64},
    B::NTuple{3,NTuple{2,Float64}},
    p::NTuple{2,Float64}
)

    alpha::Real, beta::Real = p

    z_grid, t_grid, d_grid, y_grid, P_grid, a_grid = initialise_params_test(pulse_params, N, B)

    v_grid::LinRange{Float64,Int64} = fftshift(fftfreq(length(y_grid), 1 / step(y_grid)))
    phasefactor_grid = factor_grid(z_grid, v_grid, beta)


    for i in eachindex(z_grid)
        Threads.@threads for l in eachindex(y_grid)

            evolve_polarisation_rk4!(
                view(P_grid, i, :, l),
                view(d_grid, i, :, l),
                view(a_grid, i, :, l),
                step(t_grid)
            )

            evolve_a_fft_grid_test(t_grid, a_grid, l, alpha, P_grid, step(z_grid), i, phasefactor_grid)

        end
    end

    return (P_grid, a_grid), (z_grid, t_grid, y_grid)
end