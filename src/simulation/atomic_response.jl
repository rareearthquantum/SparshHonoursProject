"""
    make_unrotate_sigma_grid(detunings, t_vec)

Return the detuning-by-time complex phase matrix `cis(-detuning * t)`
used to transform atomic coherence back before averaging.
"""
make_unrotate_sigma_grid(detunings, t_vec) = [cis(-detuning * t) for detuning in detunings, t in t_vec]

"""
    compute_polarisation_column!( polarisation_column, atomic_state, rabi_frequency_column, unrotate_sigma_grid, t_vec, rotate_sigma, cache)

Compute one transverse column of the ensemble-averaged coherence.
Reset `atomic_state[detuning, component, time]` to zero initial coherence
and inversion -1, integrate under `rabi_frequency_column`, undo detuning
rotation, and average equally over detunings into `polarisation_column`.
`unrotate_sigma_grid` has shape `(Nd, Nt)` and `rotate_sigma` contains
interleaved full- and half-time-step phases. Mutates the output column,
atomic trajectory, and RK4 cache; returns `nothing`.
"""
function compute_polarisation_column!(
    polarisation_column, atomic_state, rabi_frequency_column, unrotate_sigma_grid, t_vec, rotate_sigma, cache)
    fill!(polarisation_column, 0)

    atomic_state[:, 1, 1] .= 0.0
    atomic_state[:, 2, 1] .= -1.0
    @views rk4_custom!(atomic_state, t_vec, (rabi_frequency_column, rotate_sigma), cache)

    @inbounds for j in axes(atomic_state, 3)
        value = zero(eltype(polarisation_column))
        for d in axes(atomic_state, 1)
            value += atomic_state[d, 1, j] * unrotate_sigma_grid[d, j]
        end
        polarisation_column[j] = value / size(atomic_state, 1)
    end

    return nothing
end
