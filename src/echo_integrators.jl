mutable struct AB2Cache{T}
    fprev::T
    fcurr::T
    first_step::Bool
end

AB2Cache(u) = AB2Cache(similar(u), similar(u), true)

function ab2_step!(f, uf, ui, t, dt, p, cache::AB2Cache)
    if cache.first_step
        f(cache.fprev, ui, t, p)
        @. uf = ui + dt * cache.fprev
        cache.first_step = false
    else
        f(cache.fcurr, ui, t, p)
        @. uf = ui + (dt/2) * (3 * cache.fcurr - cache.fprev)
        copyto!(cache.fprev, cache.fcurr)
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

RK4Cache(u) = RK4Cache(similar(u), similar(u), similar(u), similar(u), similar(u))

function rk4_step!(f, uf, ui, t, dt, p, cache::RK4Cache)
    k1, k2, k3, k4 = cache.k1, cache.k2, cache.k3, cache.k4
    temp = cache.temp

    f(k1, ui, t, p)

    @. temp = ui + (dt/2) * k1
    f(k2, temp, t + dt/2, p)

    @. temp = ui + (dt/2) * k2
    f(k3, temp, t + dt/2, p)

    @. temp = ui + dt * k3
    f(k4, temp, t + dt, p)

    @. uf = ui + (dt/6) * (k1 + 2 * k2 + 2 * k3 + k4)

    return nothing
end

function rk4!(f, u, t_vec, p; substeps::Integer=1)
    substeps >= 1 || throw(ArgumentError("substeps must be at least 1"))

    dt = step(t_vec)
    cache = @views RK4Cache(u[:, 1])
    work = @views similar(u[:, 1])

    @inbounds for i in 1:(length(t_vec)-1)
        @views rk4_step!(f, u[:, i+1], u[:, i], t_vec[i], dt, p, cache)
    end

    return nothing
end


function rk4_no_substeps!(f, u, t_vec, p)
    dt = step(t_vec)
    cache = @views RK4Cache(u[:, 1])

    @inbounds for i in 1:(length(t_vec)-1)
        @views rk4_step!(f, u[:, i+1], u[:, i], t_vec[i], dt, p, cache)
    end

    return nothing
end


struct RK4whichstep
    first::Bool
    half::Bool
    last::Bool
end

function rk4_step_new!(f, uf, ui, index, dt, p, cache::RK4Cache, whichsteps)
    k1, k2, k3, k4 = cache.k1, cache.k2, cache.k3, cache.k4
    temp = cache.temp

    f(k1, ui, p, whichsteps[1], index)

    @. temp = ui + (dt/2) * k1
    f(k2, temp, p, whichsteps[2], index)

    @. temp = ui + (dt/2) * k2
    f(k3, temp, p, whichsteps[2], index)

    @. temp = ui + dt * k3
    f(k4, temp, p, whichsteps[3], index)

    @. uf = ui + (dt/6) * (k1 + 2 * k2 + 2 * k3 + k4)

    return nothing
end

function rk4_new!(f, u, t_vec, p)

    dt = step(t_vec)
    cache = @views RK4Cache(u[:, 1])
    whichsteps = (RK4whichstep(true, false, false), RK4whichstep(false, true, false), RK4whichstep(false, false, true))

    @inbounds for i in 1:(length(t_vec)-1)
        @views rk4_step_new!(f, u[:, i+1], u[:, i], i, dt, p, cache, whichsteps)
    end

    return nothing
end




#------------------------

function rk4_custom_step!(uf, ui, index, dts, p, cache::RK4Cache)
    k1, k2, k3, k4 = cache.k1, cache.k2, cache.k3, cache.k4
    temp = cache.temp

    Omega, rotate = p
    interpOmega = 0.5 * (Omega[index] + Omega[index+1])

    @views begin

        @. k1[:, 1] = -0.5im * Omega[index] * ui[:, 2] * rotate[:, 2index-1]
        @. k1[:, 2] = 2imag(conj(Omega[index]) * ui[:, 1] * conj(rotate[:, 2index-1]))

        @. temp = ui + dts[2] * k1
        @. k2[:, 1] = -0.5im * interpOmega * temp[:, 2] * rotate[:, 2index]
        @. k2[:, 2] = 2imag(conj(interpOmega) * temp[:, 1] * conj(rotate[:, 2index]))

        @. temp = ui + dts[2] * k2
        @. k3[:, 1] = -0.5im * interpOmega * temp[:, 2] * rotate[:, 2index]
        @. k3[:, 2] = 2imag(conj(interpOmega) * temp[:, 1] * conj(rotate[:, 2index]))

        @. temp = ui + dts[1] * k3
        @. k4[:, 1] = -0.5im * Omega[index+1] * temp[:, 2] * rotate[:, 2index+1]
        @. k4[:, 2] = 2imag(conj(Omega[index+1]) * temp[:, 1] * conj(rotate[:, 2index+1]))

    end

    @. uf = ui + dts[3]*(k1 + 2k2 + 2k3 + k4)

    return nothing
end


function rk4_custom!(u, time_vec, p)

    dt = step(time_vec)
    dts = (dt, dt/2, dt/6)
    cache = @views RK4Cache(u[:, :, 1])

    @inbounds for i in 1:(length(time_vec)-1)
        @views rk4_custom_step!(u[:, :, i+1], u[:, :, i], i, dts, p, cache)
    end

    return nothing
end
