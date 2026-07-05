function make_rotate_grid(detunings, time_vec)
    return [cis( δ * t) for δ in detunings, t in time_vec]
end

function initialise_field(cfg::EchoConfig, time_vec, z_vec, omega_input)
    Omega = zeros(ComplexF64, length(time_vec), length(z_vec))
    Omega[:, 1] .= omega_input.(time_vec)
    return Omega
end

function initialise_polarisation(time_vec, z_vec)
    return zeros(ComplexF64, length(time_vec), length(z_vec))
end

function reset_atom_history!(sigma_temp)
    fill!(sigma_temp, 0)
    sigma_temp[2, 1] = -1.0
    return nothing
end

function atom_stage_params(Omega_col, detuning, time_vec)
    p = (Omega_col, detuning, time_vec)
    return (p, p, p, p)
end

function atom_substeps(detunings, time_vec; max_phase_step=1.0)
    max_phase_step > 0 || throw(ArgumentError("max_phase_step must be positive"))
    maximum_phase_rate = maximum(abs, detunings)
    # Keep the fastest detuning phase advance controlled independently of Nt.
    return max(1, ceil(Int, maximum_phase_rate * step(time_vec) / max_phase_step))
end

function compute_polarisation_column!(
        P_col, sigma_temp, Omega_col, detunings, rotate_grid, time_vec, substeps)
    Nd = size(rotate_grid, 1)
    fill!(P_col, 0)

    @inbounds for i in 1:Nd
        reset_atom_history!(sigma_temp)

        rotate_vec = @view rotate_grid[i, :]
        ps = atom_stage_params(Omega_col, detunings[i], time_vec)

        rk4_staged!(atom_single!, sigma_temp, time_vec, ps; substeps=substeps)

        @views @. P_col += sigma_temp[1, :] * conj(rotate_vec)
    end

    @. P_col /= Nd

    return nothing
end

function run_propagation(cfg::EchoConfig=EchoConfig();
        omega_input=make_default_omega_input(cfg),
        compute_final_polarisation::Bool=true,
        max_atom_phase_step::Real=1.0)

    detunings = make_detunings(cfg)
    time_vec = make_time_grid(cfg)
    z_vec = make_z_grid(cfg)
    dz = step(z_vec)

    Omega = initialise_field(cfg, time_vec, z_vec, omega_input)
    P = initialise_polarisation(time_vec, z_vec)
    sigma_temp = zeros(ComplexF64, 2, length(time_vec))

    rotate_grid = make_rotate_grid(detunings, time_vec)
    substeps = atom_substeps(detunings, time_vec; max_phase_step=max_atom_phase_step)
    field_cache = @views AB2Cache(Omega[:, 1])

    @inbounds for j in 1:(length(z_vec)-1)
        @views compute_polarisation_column!(
            P[:, j], sigma_temp, Omega[:, j], detunings, rotate_grid, time_vec, substeps)

        @views ab2_step!(
            field!,
            Omega[:, j+1],
            Omega[:, j],
            z_vec[j],
            dz,
            (cfg.alpha, P[:, j]),
            field_cache,
        )
    end

    if compute_final_polarisation
        @views compute_polarisation_column!(
            P[:, end], sigma_temp, Omega[:, end], detunings, rotate_grid, time_vec, substeps)
    end

    return (
        cfg=cfg,
        detunings=detunings,
        time_vec=time_vec,
        z_vec=z_vec,
        Omega=Omega,
        P=P,
        rotate_grid=rotate_grid,
        atom_substeps=substeps,
        max_atom_phase_step=max_atom_phase_step,
    )
end
