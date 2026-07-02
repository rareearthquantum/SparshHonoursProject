using DrWatson
@quickactivate "SparshHonoursProject"

using LinearAlgebra, Plots, LinearInterpolations
include(srcdir("old/input_pulse_methods.jl"))




function atom_single!(ds, s, t, p)
    detuning, Omega, rotate_vec, index, time_vec = p
    index = time_index_grabber(t, time_vec)

    ds[1] = -0.5*im*Omega[index]*s[2]*rotate_vec[index]
    ds[2] = 2*imag(conj(Omega[index])*s[1]*conj(rotate_vec[index]))

    return nothing
end

function field!(dOmega, Omega, z, p)
    alpha, P = p

    @. dOmega = im*alpha*P
    return nothing
end

function euler_step!(f, uf, ui, t, dt, p, k)
    f(k, ui, t, p)

    @. uf = ui + dt*k

    return nothing
end

mutable struct AB2Cache{T}
    fprev::T
    fcurr::T
    first_step::Bool
end

function AB2Cache(u)
    return AB2Cache(similar(u), similar(u), true)
end

function reset_ab2_cache!(cache::AB2Cache)
    cache.first_step = true
    return nothing
end

function ab2_step!(f, uf, ui, t, dt, p, cache::AB2Cache)
    if cache.first_step
        euler_step!(f, uf, ui, t, dt, p, cache.fprev)
        cache.first_step = false
    else
        f(cache.fcurr, ui, t, p)

        @. uf = ui + (dt/2)*(3*cache.fcurr - cache.fprev)

        copyto!(cache.fprev, cache.fcurr)
    end

    return nothing
end

function ab2!(f, u, t_vec, p)
    dt = step(t_vec)
    @views cache = AB2Cache(u[:, 1])

    @inbounds for i in 1:(length(t_vec)-1)
        @views ab2_step!(f, u[:, i+1], u[:, i], t_vec[i], dt, p, cache)
    end

    return nothing
end

mutable struct RK4Cache{T}
    k1::T
    k2::T
    k3::T
    k4::T
    temp::T
end

function RK4Cache(u)
    return RK4Cache(similar(u), similar(u), similar(u), similar(u), similar(u))
end

function rk4_staged_step!(f, uf, ui, t, dt, ps, cache::RK4Cache)
    p1, p2, p3, p4 = ps
    k1, k2, k3, k4 = cache.k1, cache.k2, cache.k3, cache.k4
    temp = cache.temp

    f(k1, ui, t, p1)

    @. temp = ui + (dt/2)*k1
    f(k2, temp, t + dt/2, p2)

    @. temp = ui + (dt/2)*k2
    f(k3, temp, t + dt/2, p3)

    @. temp = ui + dt*k3
    f(k4, temp, t + dt, p4)

    @. uf = ui + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)

    return nothing
end

function rk4_staged!(f, u, t_vec, ps)
    dt = step(t_vec)
    @views cache = RK4Cache(u[:, 1])

    @inbounds for i in 1:(length(t_vec)-1)
        @views rk4_staged_step!(f, u[:, i+1], u[:, i], t_vec[i], dt, ps, cache)
    end

    return nothing
end

function time_index_grabber(t, t_vec)::Int64
    (t < t_vec[begin] || t_vec[end] < t) && return error("time out of time vec")
    return (t-Ti)÷step(t_vec)+1
end

Nd = 100
d_width = 10.0
detunings = LinRange(-d_width/2, d_width/2, Nd)

Nt = 1000
Ti = 0.0
Tf = 10.0
time_vec = LinRange(Ti, Tf, Nt)

t_width = Tf - Ti

Nz = 100
Zi = 0.0
Zf = 10.0
z_vec = LinRange(Zi, Zf, Nz)
dz = step(z_vec)

Omega = zeros(ComplexF64, Nt, Nz)
Omega_input(t) = pulse(t, 3t_width/10, t_width/10, 4pi)
Omega[:, 1] = Omega_input.(time_vec)

P = zeros(ComplexF64, Nt, Nz)
P[1, :] .= 0

sigma_temp = zeros(ComplexF64, 2, Nt)
sigma_temp[2, 1] = -1.0

field_cache = AB2Cache(Omega[:, 1])
Omega_halfsteps = similar(Omega[:, 1])

index::Int64 = 0

alpha = 1.0

rotate_grid = Array{ComplexF64}(undef, Nd, Nt)
for i in 1:Nd
    for j in 1:Nt
        rotate_grid[i, j] = exp(im*detunings[i]*time_vec[j])
    end
end

@time for j in 1:(length(z_vec)-1)

    sigma_temp[1, :] .= 0
    @. Omega_halfsteps = (Omega[:, j+1] + Omega[:, j])/2

    for i in 1:Nd
        @views rk4_staged!(atom_single!, sigma_temp, time_vec,
            (
                (detunings[i], Omega[:, j], rotate_grid[i, :], index, time_vec),
                (detunings[i], Omega[:, j], rotate_grid[i, :], index, time_vec),
                (detunings[i], Omega[:, j], rotate_grid[i, :], index, time_vec),
                (detunings[i], Omega[:, j], rotate_grid[i, :], index, time_vec)
            )
        )

        @. P[:, j] += sigma_temp[1, :]*conj(rotate_grid[i, :])
    end
    P[:, j] .*= (1/Nd)

    @views ab2_step!(
        field!,
        Omega[:, j+1],
        Omega[:, j],
        z_vec[j],
        dz,
        (alpha, P[:, j]),
        field_cache
    )
end



#PLOT
heatmap_P_abs2 = heatmap(time_vec, z_vec, abs2.(P[:, :]'), title="abs2 of polarisation", xlabel="t", ylabel="z", c=:viridis)
heatmap_P_real = heatmap(time_vec, z_vec, real.(P[:, :]'), title="real part of polarisation", xlabel="t", ylabel="z", c=:viridis)
heatmap_P_imag = heatmap(time_vec, z_vec, imag.(P[:, :]'), title="imaginary part of polarisation", xlabel="t", ylabel="z", c=:viridis)
heatmap_Omega_abs2 = heatmap(time_vec, z_vec, abs2.(Omega[:, :]'), title="abs2 of Omega", xlabel="t", ylabel="z", c=:viridis)
heatmap_Omega_real = heatmap(time_vec, z_vec, real.(Omega[:, :]'), title="real part of Omega", xlabel="t", ylabel="z", c=:viridis)
heatmap_Omega_imag = heatmap(time_vec, z_vec, imag.(Omega[:, :]'), title="imaginary part of Omega", xlabel="t", ylabel="z", c=:viridis)
plot_P_abs2 = plot(time_vec, abs2.(P[:, end÷2]), title="abs2(P) halfway through crystal")
plot_P_real = plot(time_vec, real.(P[:, end÷2]), title="real(P) halfway through crystal")
plot_P_imag = plot(time_vec, imag.(P[:, end÷2]), title="imag(P) halfway through crystal")
big_heatie = plot(heatmap_P_abs2, heatmap_P_real, heatmap_P_imag, heatmap_Omega_abs2, heatmap_Omega_real, heatmap_Omega_imag, plot_P_abs2, plot_P_real, plot_P_imag, size=(1400, 1000), layout=(3, 3))
display(big_heatie)

#=
#animation not very nice since we are still in the speed of light frame
anim = @animate for i in 1:Nt

    p1 = plot(
        z_vec,
        abs2.(Omega[i,:]),
        ylabel="Re(Ω)",
        ylim=extrema(abs2.(Omega)),
        legend=false,
        title="t = $(round(time_vec[i], digits=3))",
    )

    p2 = plot(
        z_vec,
        abs2.(P[i,:]),
        xlabel="z",
        ylabel="Re(P)",
        ylim=extrema(abs2.(P)),
        legend=false,
    )

    plot(p1, p2, layout=(2,1), size=(700,600))
end

gif(anim, "field_and_polarisation.gif", fps=30)
=#