using DrWatson
@quickactivate "SparshHonoursProject"

using LinearAlgebra
using Plots
include(srcdir("input_pulse_methods.jl"))



#=
function atom!(du, u, t, p)

    detuning, kappa, E = p

    du[1] = -detuning*u[2]
    du[2] = detuning*u[1]+kappa*E(t)*u[3]
    du[3] = -kappa*E(t)*u[2]
end
=#

function atom(u, t, p)

    detuning, kappa, E = p

    return [
        -detuning * u[2],
        detuning * u[1] + kappa * E(t) * u[3],
        -kappa * E(t) * u[2]
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
    for i in 1:(length(t_vec)-1)
        u[:, i+1] = rk4_step(f, view(u, :, i), t_vec[i], dt, p, (k1, k2, k3, k4))
    end
end

Nd = 1000
d_width = 100.0
detunings = LinRange(-d_width/2,d_width,Nd)
kappa = 1.0

Nt = 1000
Ti = 0.0
Tf = 10.0
time_vec = LinRange(Ti, Tf, Nt)

t_width = Tf - Ti


rho = Array{ComplexF64}(undef, 3, Nt, Nd)
for i in 1:Nd
    rho[:, 1, i] = [0.0, 0.0, -1.0]
end

E0(t) = 0.0
E1(t) = pulse(t, t_width/2, t_width/100, pi/2)
E2(t) = pulse(t, t_width/10, t_width/100, pi/2) + pulse(t, 4t_width/10, t_width/1000, pi/1)

for i in 1:Nd
    rk4!(atom, view(rho,:,:,i), time_vec, (detunings[i], kappa, E2))
end

sigma_minus = rho[1,:,:] + im.*rho[2,:,:]
polarisation = (1/Nd) .* vec(sum(sigma_minus,dims=2))

plot_real = plot(time_vec, abs2.(polarisation), title="abs2 of polarisation")
plot(plot_real, size=(1000, 600))