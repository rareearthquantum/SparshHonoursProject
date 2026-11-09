using FFTW

"""
    run_2d_propagation( cfg::EchoConfig=EchoConfig(); omega_2d_input=make_omega_2d_input(cfg), compute_final_polarisation::Bool=true, show_progress::Bool=true )

Solve the driven two-level ensemble with longitudinal propagation and
optional transverse diffraction. `omega_2d_input(t, y)` supplies the input
field; its default evaluates `cfg.pulses`. Set `show_progress=false` to
suppress progress output. This function does not save data or plots.

Return a named tuple containing `cfg`, `detunings`, `t_vec`, `z_vec`,
`y_vec`, `Omega[time, z, y]`, and `P[time, y]`. `P` is evaluated at the
final longitudinal plane unless `compute_final_polarisation=false`, when
it comes from the preceding plane. Local (`Nz == 1`) runs always compute
polarisation.
"""
function run_2d_propagation(
    cfg::EchoConfig=EchoConfig();
    omega_2d_input=make_omega_2d_input(cfg),
    compute_final_polarisation::Bool=true,
    show_progress::Bool=true
)
    (cfg.Nz==1) && (compute_final_polarisation=true) #can output only polarisation at Nz=1

    # Count the final atomic-response pass as work too: 100% means all requested
    # computation is finished, including local runs with no propagation steps.
    propagation_steps = cfg.Nz - 1
    total_steps = propagation_steps + Int(compute_final_polarisation)
    progress = SimulationProgress(total_steps; enabled=show_progress)

    detunings = make_detunings(cfg)
    t_vec = make_time_grid(cfg)
    z_vec = make_z_grid(cfg)
    dz = (cfg.Nz==1) ? 0.0 : step(z_vec)
    y_vec = make_y_grid(cfg)

    # Field axes: time, longitudinal position, transverse position.
    rabi_frequency = zeros(ComplexF64, length(t_vec), length(z_vec), length(y_vec))
    input_rabi_frequency = [omega_2d_input.(t_vec, y_vec[l]) for l in eachindex(y_vec)] |> stack
    rabi_frequency[:,1,:] .= input_rabi_frequency

    # Atomic-state axes: detuning, component (1=coherence, 2=inversion), time.
    atomic_state = zeros(ComplexF64, length(detunings), 2, length(t_vec))

    # Polarisation axes: time, transverse position, at the current z plane.
    polarisation = zeros(ComplexF64, length(t_vec), length(y_vec))
    polarisation_spectrum = similar(polarisation)

    unrotate_sigma_grid = make_unrotate_sigma_grid(detunings, t_vec)
    field_caches = [@views AB2Cache(rabi_frequency[:, 1, l]) for l in eachindex(y_vec)]

    ky_grid = make_ky_grid(cfg)
    diffrac_rotate = @. cis(cfg.beta*ky_grid^2*dz)
    inverse_diffrac_rotate = conj.(diffrac_rotate)

    rabi_frequency_spectrum = similar(rabi_frequency[:, 1, :])

    rotate_sigma = Array{ComplexF64}(undef, length(detunings), 2length(t_vec))
    halfdt = step(t_vec)/2
    for j in eachindex(t_vec)
        @. rotate_sigma[:, 2j-1] = cis(detunings*t_vec[j])
        @. rotate_sigma[:, 2j] = cis(detunings*(t_vec[j]+halfdt))
    end

    threadcount = Threads.nthreads()
    atom_caches = [@views RK4Cache(atomic_state[:, :, 1]) for i in 1:threadcount]

    @inbounds for j in 1:propagation_steps

        for l in eachindex(y_vec)
            @views compute_polarisation_column!(
                polarisation[:,l], atomic_state, rabi_frequency[:, j, l], unrotate_sigma_grid, t_vec, rotate_sigma, atom_caches[1])
                #need to somehow give a cache to each thread
        end

        polarisation_spectrum .= fft(polarisation, 2)
        rabi_frequency_spectrum .= fft(rabi_frequency[:, j, :], 2)

        for l in eachindex(ky_grid)
            @views ab2_step!(
                field_2d!,
                rabi_frequency[:, j+1, l],
                rabi_frequency_spectrum[:, l],
                z_vec[j],
                dz,
                (cfg.alpha, polarisation_spectrum[:, l], diffrac_rotate[l]),
                field_caches[l]
            )
        end

        for i in eachindex(t_vec)
            rabi_frequency[i, j+1, :] .*= inverse_diffrac_rotate
        end
        @views ifft!(rabi_frequency[:, j+1, :], 2)
        
        report_progress!(progress, j)

    end

    if compute_final_polarisation
        for l in eachindex(y_vec)
            # FUTURE: export requested atomic components for (last(z index), l)
            # before this column is reused, just as inside the propagation loop.
            @views compute_polarisation_column!(
                polarisation[:,l], atomic_state, rabi_frequency[:, end, l], unrotate_sigma_grid, t_vec, rotate_sigma, atom_caches[1])
        end
        # FUTURE: retain_or_write!(:P, length(z_vec), polarisation).
        report_progress!(progress, total_steps)
    end

    # Preserve the existing public result names used by plots and saved files.
    return (;
        cfg, detunings, t_vec, z_vec, y_vec,
        Omega=rabi_frequency, P=polarisation,
    )
end
