using DrWatson
@quickactivate "SparshHonoursProject"

using LinearAlgebra
using Plots
include(srcdir("input_pulse_methods.jl"))


function atom(u, t, p)

    detuning, Omega = p

    return [
        -detuning * u[2],
        detuning * u[1] + Omega(t) * u[3],
        -Omega(t) * u[2]
    ]
end

function rk4_step(f, u, t, dt, p, k)
    k1, k2, k3, k4 = k

    k1 = f(u, t, p)
    k2 = f(u + dt/2 * k1, t + dt/2, p)
    k3 = f(u + dt/2 * k2, t + dt/2, p)
    k4 = f(u + dt * k3, t + dt, p)

    return u + (dt/6) * (k1 + 2k2 + 2k3 + k4)
end

function rk4!(f, u, t_vec, p)
    k1 = similar(u[:, 1])
    k2, k3, k4 = similar(k1), similar(k1), similar(k1)
    dt = step(t_vec)
    @inbounds for i in 1:(length(t_vec)-1)
        @views u[:, i+1] = rk4_step(f, u[:, i], t_vec[i], dt, p, (k1, k2, k3, k4))
    end
end

Nd = 100
d_width = 10.0
detunings = LinRange(-d_width/2, d_width, Nd)

Nt = 1000
Ti = 0.0
Tf = 10.0
time_vec = LinRange(Ti, Tf, Nt)
t_width = Tf - Ti

rho = Array{Real}(undef, 3, Nt, Nd)
for i in 1:Nd
    rho[:, 1, i] = [0.0, 0.0, -1.0]
end

Omega(t) = pulse(t, t_width/10, t_width/100, pi/2) + pulse(t, 3t_width/10, t_width/100, pi)

@time for i in 1:Nd
    @views rk4!(atom, rho[:, :, i], time_vec, (detunings[i], Omega))
end

@show maximum(rho[1,:,:])
@show maximum(rho[2,:,:])

sigma_minus = 0.5(rho[1, :, :] - im .* rho[2, :, :]) #important s- = (1/2)(sx-isy)
polarisation = (1/Nd) .* vec(sum(sigma_minus, dims=2))
#polarisation = sigma_minus[:,1]
sum_sigma_z = (1/Nd) .* vec(sum(rho[3,:,:], dims=2))

plot_abs2 = plot(time_vec, abs2.(polarisation), title="abs2 of polarisation")
plot_real = plot(time_vec, real.(polarisation), title="real part of polarisation")
plot_imag = plot(time_vec, imag.(polarisation), title="imaginary part of polarisation")
plot_sigmaz = plot(time_vec, sum_sigma_z, title="normalised sums of sigma z")
#plot(plot_abs2, plot_real, plot_imag, ylims=(-1.1, 1.1), size=(1000, 800), layout=(3, 1))




function atom_2(u, t, p)

    detuning, Omega = p

    return [
        -im*detuning * u[1] - 0.5 * im * Omega(t) * u[2],
        im*Omega(t)*conj(u[1]) - im*conj(Omega(t))*u[1]
    ]
end

rho_2 = Array{ComplexF64}(undef, 2, Nt, Nd)
for i in 1:Nd
    rho_2[:, 1, i] = [0.0+0.0im, -1.0]
end

@time for i in 1:Nd
    @views rk4!(atom_2, rho_2[:, :, i], time_vec, (detunings[i], Omega))
end

polarisation_2 = (1/Nd) .* vec(sum(rho_2[1, :, :], dims=2))
sum_sigma_z_2 = (1/Nd) .* vec(sum(rho_2[2,:,:], dims=2))

plot_abs2_2 = plot(time_vec, abs2.(polarisation_2), title="abs2 of polarisation")
plot_real_2 = plot(time_vec, real.(polarisation_2), title="real part of polarisation")
plot_imag_2 = plot(time_vec, imag.(polarisation_2), title="imaginary part of polarisation")
plot_sigmaz_2 = plot(time_vec, real.(sum_sigma_z_2), title="normalised sums of sigma z")
plot(plot_abs2, plot_abs2_2, plot_real, plot_real_2, plot_imag, plot_imag_2, plot_sigmaz, plot_sigmaz_2, ylims=(-1.1, 1.1), size=(1000, 1000), layout=(4, 2))


#=

#=
function atom_vec(u, t, p)

    detunings, kappa, E = p
    Nd = length(detunings)

    return vcat(
        -detunings .* u[(Nd+1):2Nd],
        detunings .* u[1:Nd] .+ kappa * E(t) .* u[(2Nd+1):3Nd],
        -kappa * E(t) .* u[(Nd+1):2Nd]
    )
end
=#

function atom_vec!(du, u, t, p)

    detunings, kappa, E = p
    Nd = length(detunings)

    @. du[1:Nd] = -detunings * u[(Nd+1):2Nd]
    @. du[(Nd+1):2Nd] = detunings * u[1:Nd] + kappa * E(t) * u[(2Nd+1):3Nd]
    @. du[(2Nd+1):3Nd] = -kappa * E(t) * u[(Nd+1):2Nd]
end

function rk4_step!(f!, du, u, t, dt, p, k)
    k1, k2, k3, k4 = k

    f!(k1, u, t, p)
    f!(k2, u + dt/2 * k1, t + dt/2, p)
    f!(k3, u + dt/2 * k2, t + dt/2, p)
    f!(k4, u + dt * k3, t + dt, p)

    @. du = (dt/6) * (k1 + 2k2 + 2k3 + k4)
end

function rk4_again!(f, u, t_vec, p)
    k1 = similar(u[:, 1])
    k2, k3, k4 = similar(k1), similar(k1), similar(k1)
    dt = step(t_vec)
    @inbounds for i in 1:(length(t_vec)-1)
        @views rk4_step!(f, u[:, i+1], u[:, i], t_vec[i], dt, p, (k1, k2, k3, k4))
    end
end


rho_again = Array{ComplexF64}(undef, 3Nd, Nt)
rho_again[1:Nd, 1] .= 0.0
rho_again[(Nd+1):2Nd, 1] .= 0.0
rho_again[(2Nd+1):3Nd, 1] .= -1.0

@time rk4_again!(atom_vec!, rho_again, time_vec, (detunings, kappa, E2))

sigma_minus_again = rho_again[1:Nd, :] + im .* rho_again[(Nd+1):2Nd, :]
@show size(sigma_minus_again)
polarisation_again = (1/Nd) .* vec(sum(sigma_minus_again, dims=1))

plot_abs2_again = plot(time_vec, abs2.(polarisation_again), title="abs2 of polarisation")
plot(plot_abs2_again, size=(1000, 600))

=#