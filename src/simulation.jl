using FFTW

current_runtime(start_time) = round(time() - start_time; digits=2)

make_unrotate_sigma_grid(detunings, t_vec) = [cis(-detuning * t) for detuning in detunings, t in t_vec]

function compute_polarisation_column!(
    polarisation_column, atomic_state, rabi_frequency_column, unrotate_sigma_grid, t_vec, rotate_sigma, cache)
    # STORAGE DESIGN ONLY: a future atomic_output callback would receive coherence
    # and inversion for this (z_index, y_index) before atomic_state is reused.
    # The caller must supply both indices; this buffer holds ONE transverse column,
    # not the atomic state of the whole transverse plane.
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

function run_2d_propagation(
    cfg::EchoConfig=EchoConfig();
    omega_2d_input=make_omega_2d_input(cfg),
    compute_final_polarisation::Bool=true,
    show_progress::Bool=true
)
    # STORAGE DESIGN ONLY -- the following pseudocode is not implemented.
    # Proposed extra keyword, independent of the numerical configuration:
    # storage = (ram_fields=(:Omega,), disk_fields=(), z_stride=1,
    #            directory=nothing, checkpoint_stride=0, buffer_planes=1)
    # Available history fields: :Omega, :P, :coherence, :inversion. Each may be
    # retained in RAM, written to disk, both, or neither. Required working buffers
    # are still computed regardless of these choices. ram_fields=() means a rolling
    # workspace; it does not mean that the calculation uses no RAM.
    # The default above describes today's Omega history plus final P result.
    # Validate names, positive strides/buffer sizes, and a writable output path
    # when disk output or checkpoints are requested. Checkpoint stride 0 disables
    # restart files. Estimate RAM and disk requirements before allocating: one
    # complex field plane costs 16*Nt*Ny bytes, an atomic component history can cost
    # 16*Nd*Nt*Ny*saved_z_count bytes. Compression should be optional and lossless.
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

    # FUTURE RAM/DISK SETUP: replace the full Omega allocation below with two
    # working planes and histories allocated only for storage.ram_fields.
    # workspace = (current=zeros(ComplexF64, cfg.Nt, cfg.Ny),
    #              next=zeros(ComplexF64, cfg.Nt, cfg.Ny))
    # workspace.current .= omega_2d_input.(t_vec, transpose(y_vec))
    # history = allocate_selected_histories(storage.ram_fields, cfg)
    # archive = needs_disk_output(storage) ?
    #     open_new_run(storage.directory; cfg, detunings, t_vec, z_vec, y_vec) : nothing
    # open_new_run would also record schema version, axes, numeric types, unit
    # conventions, code/environment versions, storage policy, and input identity.
    # A custom omega_2d_input needs a retained recipe/source too; cfg alone does
    # not describe an arbitrary callable well enough to reproduce or resume it.
    # Use a new run directory; never overwrite an existing run implicitly.
    # try ... finally (isnothing(archive) || close(archive)) end must cover the run.
    # Pseudocode helpers in these comments are future interfaces, not existing APIs.
    # Keep full (t,y) working planes: atomic integration needs the time trace and
    # the transverse FFT needs all y samples. Disk output is organised by z plane.

    # Field axes: time, longitudinal position, transverse position.
    rabi_frequency = zeros(ComplexF64, length(t_vec), length(z_vec), length(y_vec))
    input_rabi_frequency = [omega_2d_input.(t_vec, y_vec[l]) for l in eachindex(y_vec)] |> stack
    rabi_frequency[:,1,:] .= input_rabi_frequency

    # Atomic-state axes: detuning, component (1=coherence, 2=inversion), time.
    # FUTURE ATOMIC RAM POLICY: this Nd x 2 x Nt workspace is already reused across
    # y and z; do not allocate an Nd x 2 x Nt x Ny x Nz atomic history by default.
    # A further refactor of rk4_custom! could retain only current/next Nd x 2
    # states, accumulating polarisation and streaming selected atomic time blocks.
    # That needs a new RK4/output interface; simply shrinking this array would
    # break compute_polarisation_column! as currently written. The precomputed
    # rotation grids also remain O(Nd*Nt) unless phases are generated as needed.
    atomic_state = zeros(ComplexF64, length(detunings), 2, length(t_vec))

    # Polarisation axes: time, transverse position, at the current z plane.
    polarisation = zeros(ComplexF64, length(t_vec), length(y_vec))
    polarisation_spectrum = similar(polarisation)

    unrotate_sigma_grid = make_unrotate_sigma_grid(detunings, t_vec)
    field_caches = [@views AB2Cache(rabi_frequency[:, 1, l]) for l in eachindex(y_vec)]
    # FUTURE ROLLING MODE: construct these AB2 caches from workspace.current.
    # Keep their previous derivative and first_step flag across ALL z steps;
    # resetting them when rotating buffers would change the integration method.

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

    # FUTURE SHARED OUTPUT HOOK: retain_or_write!(:Omega, 1, workspace.current).
    # The hook applies ram_fields and disk_fields independently. Store the input
    # and final planes for selected fields, when computed, even when z_stride
    # skips intermediate output planes.
    # z_stride changes output sampling only; every propagation step still runs.
    # A simple disk layout is run/Omega/z000001.jld2, run/P/z000001.jld2, etc.,
    # containing one plane and its true z index. Atomic output additionally needs
    # y indices and detuning/time axes, with bounded tiles if a column is too large.
    # Stage each write, then publish a completion marker after a successful write;
    # readers should ignore incomplete files. Avoid rewriting a full run file on
    # every step, or repeatedly jldsave-ing different planes to the same filename.

    @inbounds for j in 1:propagation_steps

        for l in eachindex(y_vec)
            # FUTURE: read the field from workspace.current[:, l]. After RK4,
            # retain_or_write_atomic!(j, l, atomic_state) must run here, BEFORE the
            # next column overwrites it. Persist only the requested components;
            # coherence is stored in the rotating frame used by atomic_state.
            # Record that convention, or explicitly unrotate it for exported data.
            @views compute_polarisation_column!(
                polarisation[:,l], atomic_state, rabi_frequency[:, j, l], unrotate_sigma_grid, t_vec, rotate_sigma, atom_caches[1])
                #need to somehow give a cache to each thread
        end

        # FUTURE: retain_or_write!(:P, j, polarisation) BEFORE reusing P.
        # This is polarisation at z[j], not at the new field plane z[j+1].
        # An optional diagnostics hook can accumulate energy, area, and selected
        # traces here, allowing small summary outputs without full field histories.
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

        # FUTURE: AB2 writes workspace.next in Fourier space; apply the diffraction
        # factor and inverse FFT exactly as above before exporting the physical
        # field: retain_or_write!(:Omega, j + 1, workspace.next).
        # Then swap: workspace = (current=workspace.next, next=workspace.current).
        # No earlier field plane needs to remain in RAM. Keep P, FFT workspaces,
        # atomic workspaces, and integrator caches; they are required for stepping.
        # Synchronous writes are simplest. Any asynchronous writer needs a bounded
        # queue and an owned copy: views into these buffers become invalid on reuse.
        # Propagate writer failures, apply backpressure, and drain writes on exit.

        # FUTURE CHECKPOINT: after a completed step and buffer swap, optionally save
        # (field_z_index=j+1, current_field, AB2 fprev/first_step per ky, cfg/grids,
        #  archive progress per field, input/code identity, elapsed time).
        # P is still at index j here; record that independently. Output strides
        # must never determine restart position. On resume, verify configuration
        # and input compatibility, rebuild scratch/phase/RK4 buffers, restore AB2
        # history, and continue at the next step without duplicating saved planes.
        # Write checkpoints separately from sparse scientific output, publishing
        # only a complete checkpoint after pending output has been committed.
        
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

    # FUTURE FINALISATION: Nz=1 never enters the propagation loop; still export its
    # input/final field once and its local P/atomic output via the final calculation.
    # If compute_final_polarisation=false, P belongs to the preceding z plane;
    # record that index and the missing final-P sample rather than relabelling it.
    # Drain output, mark the archive complete, and close it in the finally block.
    # On failure, keep committed planes/checkpoints and mark the run incomplete.
    # Rolling-mode results should contain final field snapshots, retained histories
    # (if any), and archive references with saved indices. Copy the final workspace
    # before handing it to a caller if another task may still reuse that buffer.
    # Existing plotting/output code expects result.Omega[t,z,y] in RAM, so provide
    # explicit read_field(:Omega; z_indices=...) / read_field(:P; z_indices=...)
    # adapters and slab-wise plotting/reductions. Never silently reload the entire
    # archive, substitute zeros for unsaved fields, or duplicate disk-only data via
    # save_result. Missing histories should produce a clear request for the needed
    # saved planes. Check saved indices/cadence before plotting or curve fitting.

    # Preserve the existing public result names used by plots and saved files.
    return (;
        cfg, detunings, t_vec, z_vec, y_vec,
        Omega=rabi_frequency, P=polarisation,
    )
end

