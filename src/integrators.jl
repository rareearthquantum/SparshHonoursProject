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



#------------------------

function rk4_custom_step!(uf, ui, index, dts, p, cache::RK4Cache)
    k1, k2, k3, k4 = cache.k1, cache.k2, cache.k3, cache.k4
    temp = cache.temp

    rabi_frequency, rotate = p
    midpoint_rabi_frequency = 0.5 * (rabi_frequency[index] + rabi_frequency[index+1])

    @views begin

        @. k1[:, 1] = -0.5im * rabi_frequency[index] * ui[:, 2] * rotate[:, 2index-1]
        @. k1[:, 2] = 2imag(conj(rabi_frequency[index]) * ui[:, 1] * conj(rotate[:, 2index-1]))

        @. temp = ui + dts[2] * k1
        @. k2[:, 1] = -0.5im * midpoint_rabi_frequency * temp[:, 2] * rotate[:, 2index]
        @. k2[:, 2] = 2imag(conj(midpoint_rabi_frequency) * temp[:, 1] * conj(rotate[:, 2index]))

        @. temp = ui + dts[2] * k2
        @. k3[:, 1] = -0.5im * midpoint_rabi_frequency * temp[:, 2] * rotate[:, 2index]
        @. k3[:, 2] = 2imag(conj(midpoint_rabi_frequency) * temp[:, 1] * conj(rotate[:, 2index]))

        @. temp = ui + dts[1] * k3
        @. k4[:, 1] = -0.5im * rabi_frequency[index+1] * temp[:, 2] * rotate[:, 2index+1]
        @. k4[:, 2] = 2imag(conj(rabi_frequency[index+1]) * temp[:, 1] * conj(rotate[:, 2index+1]))

    end

    @. uf = ui + dts[3]*(k1 + 2k2 + 2k3 + k4)

    return nothing
end


function rk4_custom!(u, t_vec, p, cache)

    dt = step(t_vec)
    dts = (dt, dt/2, dt/6)

    @inbounds for i in 1:(length(t_vec)-1)
        @views rk4_custom_step!(u[:, :, i+1], u[:, :, i], i, dts, p, cache)
    end

    return nothing
end


function field_2d!(rabi_frequency_derivative, rabi_frequency_spectrum, z, p)
    alpha, polarisation_spectrum, rotfactor = p

    @. rabi_frequency_derivative = im * alpha * polarisation_spectrum * rotfactor

    return nothing
end