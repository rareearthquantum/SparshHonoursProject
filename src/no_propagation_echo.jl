import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using LinearAlgebra
using Plots
include(joinpath(@__DIR__, "input_pulse_methods.jl"))

#FUNCTIONS
function atom!(ds, s, t, p)

    detunings, Nd, Omega = p

    @. ds[1:Nd] = -0.5*im*Omega(t)*s[(Nd+1):2Nd]*cis(detunings*t)
    @. ds[(Nd+1):2Nd] = 2*imag(conj(Omega(t))*s[1:Nd]*cis(-detunings*t))

    return nothing
end

function rk4_step!(f, uf, ui, t, dt, p, k)
    k1, k2, k3, k4 = k

    f(k1, ui, t, p)
    f(k2, ui + dt/2 * k1, t + dt/2, p)
    f(k3, ui + dt/2 * k2, t + dt/2, p)
    f(k4, ui + dt * k3, t + dt, p)

    @. uf = ui + (dt/6) * (k1 + 2k2 + 2k3 + k4)

    return nothing
end

function rk4!(f, u, t_vec, p)
    k1 = similar(u[:, 1])
    k2, k3, k4 = similar(k1), similar(k1), similar(k1)
    dt = step(t_vec)
    @inbounds for i in 1:(length(t_vec)-1)
        @views rk4_step!(f, u[:, i+1], u[:, i], t_vec[i], dt, p, (k1, k2, k3, k4))
    end

    return nothing
end

#PARAMS
Nt = 1000
Ti = 0.0
Tf = 10.0
time_vec = LinRange(Ti, Tf, Nt)

d_width = 1/step(time_vec)
Nd = 2*d_width |> x->round(Int,x) 
@show d_width
detunings = LinRange(-d_width/2, d_width/2, Nd)

#INTIALISE
rho = Array{ComplexF64}(undef, 2Nd, Nt)
rho[1:Nd,1] .= 0.0 + 0.0im
rho[Nd+1:2Nd,1] .= -1.0

#PULSE WITH PARAMS
t_width = Tf - Ti
Omega_input(t) = pulse(t, t_width/10, t_width/100, pi/2) + pulse(t, 3t_width/10, t_width/100, pi)


#SOLVE
@time @views rk4!(atom!, rho, time_vec, (detunings, length(detunings), Omega_input))


#unrotate
for i in 1:Nt
    @. rho[1:Nd, i] *= cis(-detunings*time_vec[i])
end

#SUM AND INTERPRET
polarisation = (1/Nd) .* vec(sum(rho[1:Nd, :], dims=1))
sum_sigma_z = (1/Nd) .* vec(sum(rho[(Nd+1):2Nd, :], dims=1))

#PLOT
plot_abs2 = plot(time_vec, abs2.(polarisation), title="abs2 of polarisation")
plot_real = plot(time_vec, real.(polarisation), title="real part of polarisation")
plot_imag = plot(time_vec, imag.(polarisation), title="imaginary part of polarisation")
plot_sigmaz = plot(time_vec, real.(sum_sigma_z), title="normalised sums of sigma z")
plot(plot_abs2, plot_real, plot_imag, plot_sigmaz, size=(500, 800), layout=(4, 1))
