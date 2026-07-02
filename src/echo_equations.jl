# Right-hand sides for the atom and field equations.

function time_index_grabber(t, t_vec)::Int
    # For an array-valued field Ω(tᵢ), use the most recent grid value.
    # This preserves the behaviour of the original floor-division version,
    # but avoids depending on a global Ti variable.
    T = float(eltype(t_vec))
    tol = 100 * eps(T) * max(abs(first(t_vec)), abs(last(t_vec)), one(T))

    if t < first(t_vec) - tol || t > last(t_vec) + tol
        error("time $t is outside the time grid [$(first(t_vec)), $(last(t_vec))]")
    end

    t_clamped = clamp(t, first(t_vec), last(t_vec))
    return searchsortedlast(t_vec, t_clamped)
end

function atom_single!(ds, s, t, p)
    Omega, rotate_vec, time_vec = p
    index = time_index_grabber(t, time_vec)

    ds[1] = -0.5im * Omega[index] * s[2] * rotate_vec[index]
    ds[2] = 2imag(conj(Omega[index]) * s[1] * conj(rotate_vec[index]))

    return nothing
end

function field!(dOmega, Omega, z, p)
    alpha, P = p

    @. dOmega = im * alpha * P

    return nothing
end
