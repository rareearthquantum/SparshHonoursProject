mutable struct RK4Cache{T}
    k1::T
    k2::T
    k3::T
    k4::T
    temp::T
end

"""
    RK4Cache(u)

Allocate four stage buffers and one temporary buffer matching the shape
and element type of `u`, for reuse by the atomic RK4 integrator.
"""
RK4Cache(u) = RK4Cache(similar(u), similar(u), similar(u), similar(u), similar(u))



#------------------------

"""
    rk4_custom_step!(uf, ui, index, dts, p, cache::RK4Cache)

Advance a detuning-by-component atomic state from `ui` into `uf` using
RK4. Component 1 is coherence and component 2 is inversion.
`dts = (dt, dt/2, dt/6)` and `p = (rabi_frequency, rotate)`, where the field
is sampled on the time grid and `rotate` on its interleaved half steps.
`index` selects the starting time sample. Mutates `uf` and `cache`;
returns `nothing`.
"""
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

"""
    rk4_custom!(u, t_vec, p, cache)

Fill the atomic trajectory `u[detuning, component, time]` after its first
time slice, using the uniform range `t_vec` and reusable RK4 `cache`.
`p` contains the time-sampled field and half-step rotation grid as in
[`rk4_custom_step!`](@ref). Returns `nothing`.
"""
function rk4_custom!(u, t_vec, p, cache)

    dt = step(t_vec)
    dts = (dt, dt/2, dt/6)

    @inbounds for i in 1:(length(t_vec)-1)
        @views rk4_custom_step!(u[:, :, i+1], u[:, :, i], i, dts, p, cache)
    end

    return nothing
end
