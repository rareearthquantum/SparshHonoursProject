using DrWatson
@quickactivate "SparshHonoursProject"

using LinearAlgebra
using Plots
include(srcdir("input_pulse_methods.jl"))

#FUNCTIONS
function atom_rotate(u, t, p)

    detuning, Omega = p

    return [
        -0.5*im*Omega(t)*u[2]*exp(im*detuning*t),
        im*Omega(t)*conj(u[1]) - im*conj(Omega(t))*u[1]
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

#PARAMS
Nd = 100
d_width = 1.0
detunings = LinRange(-d_width/2, d_width, Nd)

Nt = 1000
Ti = 0.0
Tf = 10.0
time_vec = LinRange(Ti, Tf, Nt)

t_width = Tf - Ti

#INTIALISE
rho = Array{ComplexF64}(undef, 2, Nt, Nd)
for i in 1:Nd
    rho[:, 1, i] = [0.0+0.0im, -1.0]
end

#PULSE WITH PARAMS
Omega(t) = pulse(t, t_width/10, t_width/100, pi/2) + pulse(t, 3t_width/10, t_width/100, pi)


#SOLVE
@time Threads.@threads for i in 1:Nd
    @views rk4!(atom_rotate, rho[:, :, i], time_vec, (detunings[i], Omega))
end

#unrotate
for i in 1:Nt
    for j in 1:Nd
        rho[1, i, j] = rho[1, i, j]*exp(-im*detunings[j]*time_vec[i])
    end
end

#SUM AND INTERPRET
polarisation = (1/Nd) .* vec(sum(rho[1, :, :], dims=2))
sum_sigma_z = (1/Nd) .* vec(sum(rho[2,:,:], dims=2))

#PLOT
plot_abs2 = plot(time_vec, abs2.(polarisation), title="abs2 of polarisation")
plot_real = plot(time_vec, real.(polarisation), title="real part of polarisation")
plot_imag = plot(time_vec, imag.(polarisation), title="imaginary part of polarisation")
plot_sigmaz = plot(time_vec, real.(sum_sigma_z), title="normalised sums of sigma z")
plot(plot_abs2, plot_real, plot_imag, plot_sigmaz, ylims=(-1.1,1.1), size=(500, 800), layout=(4, 1))