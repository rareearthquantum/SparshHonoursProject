make_unrotate_grid(detunings, time_vec) = [cis(-detuning * t) for detuning in detunings, t in time_vec]

function compute_polarisation_column!(
        P_col, sigma_temp, Omega_col, detunings, unrotate_grid, time_vec, substeps)
    fill!(P_col, 0)

    @inbounds for i in eachindex(detunings)
        fill!(sigma_temp, 0)
        sigma_temp[2, 1] = -1.0
        rk4!(atom_single!, sigma_temp, time_vec,
            (Omega_col, detunings[i], time_vec); substeps=substeps)
        @views @. P_col += sigma_temp[1, :] * unrotate_grid[i, :]
    end

    P_col ./= length(detunings)

    return nothing
end

function run_propagation(cfg::EchoConfig=EchoConfig();
        omega_input=make_omega_input(cfg),
        compute_final_polarisation::Bool=true,
        max_atom_phase_step::Real=1.0)

    detunings = make_detunings(cfg)
    time_vec = make_time_grid(cfg)
    z_vec = make_z_grid(cfg)
    dz = step(z_vec)

    Omega = zeros(ComplexF64, length(time_vec), length(z_vec), length(y_vec))
    Omega[:, 1] .= omega_input.(time_vec)
    P = zeros(ComplexF64, length(time_vec), length(z_vec))
    sigma_temp = zeros(ComplexF64, 2, length(time_vec))

    unrotate_grid = make_unrotate_grid(detunings, time_vec)
    substeps = atom_substeps(detunings, time_vec; max_phase_step=max_atom_phase_step)
    field_cache = @views AB2Cache(Omega[:, 1])

    @inbounds for j in 1:(length(z_vec)-1)
        @views compute_polarisation_column!(
            P[:, j], sigma_temp, Omega[:, j], detunings, unrotate_grid, time_vec, substeps)

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
            P[:, end], sigma_temp, Omega[:, end], detunings, unrotate_grid, time_vec, substeps)
    end

    return (;
        cfg, detunings, time_vec, z_vec, Omega, P,
        atom_substeps=substeps,
        max_atom_phase_step,
    )
end
