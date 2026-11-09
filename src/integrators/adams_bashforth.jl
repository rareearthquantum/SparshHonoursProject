mutable struct AB2Cache{T}
    fprev::T
    fcurr::T
    first_step::Bool
end

"""
    AB2Cache(u)

Allocate derivative buffers matching `u` and mark the next AB2 step as
the first step. Use a separate cache for each independent trajectory.
"""
AB2Cache(u) = AB2Cache(similar(u), similar(u), true)

"""
    ab2_step!(f, uf, ui, t, dt, p, cache::AB2Cache)

Advance `ui` by `dt` into `uf` using two-step Adams–Bashforth, with an
Euler start on the first call. `f(du, u, t, p)` must fill `du` in place.
Mutates `uf` and `cache`, assumes a fixed step size, and returns `nothing`.
"""
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
