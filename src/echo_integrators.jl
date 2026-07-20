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


mutable struct RK4whichstep
    first::Bool
    half::Bool
    last::Bool
end

function rk4_step_new!(f, uf, ui, index, dt, p, cache::RK4Cache)
    k1, k2, k3, k4 = cache.k1, cache.k2, cache.k3, cache.k4
    temp = cache.temp

    f(k1, ui, p, RK4whichstep(true,false,false), index)

    @. temp = ui + (dt/2) * k1
    f(k2, temp, p, RK4whichstep(false,true,false), index)

    @. temp = ui + (dt/2) * k2
    f(k3, temp, p, RK4whichstep(false,true,false), index)

    @. temp = ui + dt * k3
    f(k4, temp, p, RK4whichstep(false,false,true), index)

    @. uf = ui + (dt/6) * (k1 + 2 * k2 + 2 * k3 + k4)

    return nothing
end

function rk4_new!(f, u, t_vec, p)

    dt = step(t_vec)
    cache = @views RK4Cache(u[:, 1])

    @inbounds for i in 1:(length(t_vec)-1)
        @views rk4_step_new!(f, u[:, i+1], u[:, i], i, dt, p, cache)
    end

    return nothing
end