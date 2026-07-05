function euler_step!(f, uf, ui, t, dt, p, k)
    f(k, ui, t, p)
    @. uf = ui + dt * k
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
        @. uf = ui + (dt/2) * (3 * cache.fcurr - cache.fprev)
        copyto!(cache.fprev, cache.fcurr)
    end

    return nothing
end

function ab2!(f, u, t_vec, p)
    dt = step(t_vec)
    cache = @views AB2Cache(u[:, 1])

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

    @. temp = ui + (dt/2) * k1
    f(k2, temp, t + dt/2, p2)

    @. temp = ui + (dt/2) * k2
    f(k3, temp, t + dt/2, p3)

    @. temp = ui + dt * k3
    f(k4, temp, t + dt, p4)

    @. uf = ui + (dt/6) * (k1 + 2 * k2 + 2 * k3 + k4)

    return nothing
end

function rk4_staged!(f, u, t_vec, ps; substeps::Integer=1)
    substeps >= 1 || throw(ArgumentError("substeps must be at least 1"))

    dt = step(t_vec) / substeps
    cache = @views RK4Cache(u[:, 1])
    work = @views similar(u[:, 1])

    @inbounds for i in 1:(length(t_vec)-1)
        @views copyto!(work, u[:, i])

        for substep in 1:(substeps-1)
            t = t_vec[i] + (substep - 1) * dt
            rk4_staged_step!(f, work, work, t, dt, ps, cache)
        end

        t = t_vec[i] + (substeps - 1) * dt
        @views rk4_staged_step!(f, u[:, i+1], work, t, dt, ps, cache)
    end

    return nothing
end
