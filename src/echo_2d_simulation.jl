using FFTW

make_unrotate_grid(detunings, time_vec) = [cis(-detuning * t) for detuning in detunings, t in time_vec]

function atom_substeps(detunings, time_vec; max_phase_step=1.0)
    max_phase_step > 0 || throw(ArgumentError("max_phase_step must be positive"))
    # Keep the fastest detuning phase advance controlled independently of Nt.
    return max(1, ceil(Int, maximum(abs, detunings) * step(time_vec) / max_phase_step))
end

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

function run_2d_propagation(cfg::EchoConfig=EchoConfig();
    omega_2d_input=make_omega_2d_input(cfg),
    compute_final_polarisation::Bool=true,
    max_atom_phase_step::Real=1.0)

    detunings = make_detunings(cfg)
    time_vec = make_time_grid(cfg)
    z_vec = make_z_grid(cfg)
    dz = step(z_vec)
    y_vec = make_y_grid(cfg)

    Omega = zeros(ComplexF64, length(time_vec), length(z_vec), length(y_vec))
    for l in eachindex(y_vec)
        Omega[:, 1, l] .= omega_2d_input.(time_vec, y_vec[l])
    end
    Omega[:,1,:] .= fftshift(fft(Omega[:,1,:],2),2)

    P = zeros(ComplexF64, length(time_vec), length(z_vec), length(y_vec))
    sigma_temp = zeros(ComplexF64, 2, length(time_vec))

    unrotate_grid = make_unrotate_grid(detunings, time_vec)
    substeps = atom_substeps(detunings, time_vec; max_phase_step=max_atom_phase_step)
    field_caches = [@views AB2Cache(Omega[:, 1, l]) for l in eachindex(y_vec)]

    ky_grid = 2pi .* fftshift(fftfreq(length(y_vec), 1/step(y_vec)))
    rotfactorgrid = [cis(cfg.beta*ky^2*z) for z in z_vec, ky in ky_grid]

    P_ky = similar(P[:,1,:])

    @inbounds for j in 1:(length(z_vec)-1)
        for l in eachindex(y_vec)
            @views compute_polarisation_column!(
                P[:, j, l], sigma_temp, Omega[:, j, l], detunings, unrotate_grid, time_vec, substeps)
        end

        P_ky[:,:] = fftshift(fft(P[:,j,:],2),2)

        for l in eachindex(ky_grid)
            @views ab2_step!(
                field_2d!,
                Omega[:, j+1, l],
                Omega[:, j, l],
                z_vec[j],
                dz,
                (cfg.alpha, P_ky[:, l], rotfactorgrid[j,l]),
                field_caches[l],
            )
        end
    end

    if compute_final_polarisation
        for l in eachindex(y_vec)
            @views compute_polarisation_column!(
                P[:, end, l], sigma_temp, Omega[:, end, l], detunings, unrotate_grid, time_vec, substeps)
        end
    end

    inverse_rotfactorgrid = conj.(rotfactorgrid)
    for j in eachindex(z_vec)
        for i in eachindex(time_vec)
            @views Omega[i, j, :] .*= inverse_rotfactorgrid[j, :]
        end
    end
    ifft!(Omega, 3)

    return (;
        cfg, detunings, time_vec, z_vec, y_vec, Omega, P,
        atom_substeps=substeps,
        max_atom_phase_step,
    )
end
