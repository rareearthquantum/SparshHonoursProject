function time_index_grabber(t, t_vec)::Int
    T = float(eltype(t_vec))
    tol = 100 * eps(T) * max(abs(first(t_vec)), abs(last(t_vec)), one(T))

    if t < first(t_vec) - tol || t > last(t_vec) + tol
        error("time $t is outside the time grid [$(first(t_vec)), $(last(t_vec))]")
    end

    t_clamped = clamp(t, first(t_vec), last(t_vec))
    return searchsortedlast(t_vec, t_clamped)
end

function linear_sample(values, t, time_vec)
    index = time_index_grabber(t, time_vec)
    index == lastindex(time_vec) && return values[index]

    fraction = (t - time_vec[index]) / (time_vec[index+1] - time_vec[index])
    return (1 - fraction) * values[index] + fraction * values[index+1]
end

function atom_single!(ds, s, t, p)
    Omega, detuning, time_vec = p
    # RK4 stages generally fall between saved time samples.
    Omega_t = linear_sample(Omega, t, time_vec)
    rotate_t = cis(detuning * t)

    ds[1] = -0.5im * Omega_t * s[2] * rotate_t
    ds[2] = 2imag(conj(Omega_t) * s[1] * conj(rotate_t))

    return nothing
end

function field!(dOmega, Omega, z, p)
    alpha, P = p

    @. dOmega = im * alpha * P

    return nothing
end



function interp_half(values, vec, index)
    index == lastindex(vec) && return values[index]

    return 0.5 * values[index] + 0.5 * values[index+1]
end


function field_2d!(dOmega, Omega, z, p)
    alpha, P_ky, rotfactor = p

    @. dOmega = im * alpha * P_ky * rotfactor

    return nothing
end

function atom_woah!(ds, s, t, p)
    Omega, time_vec, rotate = p

    index = searchsortedlast(time_vec, t)
    if (index == 0)
        index = searchsortedlast(time_vec, t - step(time_vec)/2)
        Omega_t = interp_half(Omega, time_vec, index)
        rotate_t = rotate[2index]
    else
        Omega_t = Omega[index]
        rotate_t = rotate[2index-1]
    end

    ds[1] = -0.5im * Omega_t * s[2] * rotate_t
    ds[2] = 2imag(conj(Omega_t) * s[1] * conj(rotate_t))

    return nothing
end