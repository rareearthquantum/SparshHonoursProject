using FFTW

current_runtime(start_time) = round(time() - start_time; digits=2)

make_unrotate_sigma_grid(detunings, time_vec) = [cis(-detuning * t) for detuning in detunings, t in time_vec]

function compute_polarisation_column!(
    P_col, sigma_temp, Omega_col, unrotate_sigma_grid, time_vec, rotate_sigma, cache)
    fill!(P_col, 0)

    sigma_temp[:, 1, 1] .= 0.0
    sigma_temp[:, 2, 1] .= -1.0
    @views rk4_custom!(sigma_temp, time_vec, (Omega_col, rotate_sigma), cache)

    @inbounds for j in axes(sigma_temp, 3)
        value = zero(eltype(P_col))
        for d in axes(sigma_temp, 1)
            value += sigma_temp[d, 1, j] * unrotate_sigma_grid[d, j]
        end
        P_col[j] = value / size(sigma_temp, 1)
    end

    return nothing
end

function run_2d_propagation(
    cfg::EchoConfig=EchoConfig();
    omega_2d_input=make_omega_2d_input(cfg),
    compute_final_polarisation::Bool=true
)
    (cfg.Nz==1) && (compute_final_polarisation=true) #can output only polarisation at Nz=1

    start_time = time()

    detunings = make_detunings(cfg)
    time_vec = make_time_grid(cfg)
    z_vec = make_z_grid(cfg)
    dz = (cfg.Nz==1) ? 0.0 : step(z_vec)
    y_vec = make_y_grid(cfg)

    Omega = zeros(ComplexF64, length(time_vec), length(z_vec), length(y_vec))
    Omega_input = [omega_2d_input.(time_vec, y_vec[l]) for l in eachindex(y_vec)] |> stack
    Omega[:,1,:] .= Omega_input

    sigma_temp = zeros(ComplexF64, length(detunings), 2, length(time_vec))

    P = zeros(ComplexF64, length(time_vec), length(y_vec))
    P_ky = similar(P)

    unrotate_sigma_grid = make_unrotate_sigma_grid(detunings, time_vec)
    field_caches = [@views AB2Cache(Omega[:, 1, l]) for l in eachindex(y_vec)]

    ky_grid = make_ky_grid(cfg)
    diffrac_rotate = @. cis(cfg.beta*ky_grid^2*dz)
    inverse_diffrac_rotate = conj.(diffrac_rotate)

    Omega_ky = similar(Omega[:, 1, :])

    rotate_sigma = Array{ComplexF64}(undef, length(detunings), 2length(time_vec))
    halfdt = step(time_vec)/2
    for j in eachindex(time_vec)
        @. rotate_sigma[:, 2j-1] = cis(detunings*time_vec[j])
        @. rotate_sigma[:, 2j] = cis(detunings*(time_vec[j]+halfdt))
    end

    threadcount = Threads.nthreads()
    atom_caches = [@views RK4Cache(sigma_temp[:, :, 1]) for i in 1:threadcount]

    percent_count = 0
    tenthofNz = max(1,length(z_vec)÷10)

    @inbounds for j in 1:(length(z_vec)-1)

        for l in eachindex(y_vec)
            @views compute_polarisation_column!(
                P[:,l], sigma_temp, Omega[:, j, l], unrotate_sigma_grid, time_vec, rotate_sigma, atom_caches[1]) 
                #need to somehow give a cache to each thread
        end

        P_ky .= fft(P, 2)
        Omega_ky .= fft(Omega[:, j, :], 2)

        for l in eachindex(ky_grid)
            @views ab2_step!(
                field_2d!,
                Omega[:, j+1, l],
                Omega_ky[:, l],
                z_vec[j],
                dz,
                (cfg.alpha, P_ky[:, l], diffrac_rotate[l]),
                field_caches[l]
            )
        end

        for i in eachindex(time_vec)
            Omega[i, j+1, :] .*= inverse_diffrac_rotate
        end
        @views ifft!(Omega[:, j+1, :], 2)

        
        if (j%tenthofNz==0)
            percent_count += 10
            println("$(percent_count)% complete in $(current_runtime(start_time))s")
        end

    end

    if compute_final_polarisation
        for l in eachindex(y_vec)
            @views compute_polarisation_column!(
                P[:,l], sigma_temp, Omega[:, end, l], unrotate_sigma_grid, time_vec, rotate_sigma, atom_caches[1])
        end
    end

    return (;
        cfg, detunings, time_vec, z_vec, y_vec, Omega, P,
    )
end


function run_simulation(cfg)
    # RUN
    @time result = run_2d_propagation(cfg)

    # Setting info for saving into file
    timestamp = Dates.format(now(), dateformat"yyyymmdd-HHMMSS-sss")
    parameter_info = "_Nt=$(cfg.Nt)_Nd=$(cfg.Nd)_Nz=$(cfg.Nz)_Ny=$(cfg.Ny)_dwidth=$(cfg.d_width)_alpha=$(cfg.alpha)_beta=$(cfg.beta)_pulsecount=$(length(cfg.pulses))_"
    for pulse in cfg.pulses
        parameter_info *= "_area=$(pulse[1].area)_width=$(pulse[1].width)_"
    end


    # Setting plot directory and helper function
    plot_output_dir = joinpath(dirname(@__DIR__), "plots", "prop_2d")
    mkpath(plot_output_dir)
    plot_n_save(func, name) = save_plot(result, func, plot_output_dir, name; parameter_info=parameter_info, timestamp=timestamp)


    # Plotting
    #plot_n_save(plot_sum_omega, "energy")
    #plot_n_save(x -> plot_sum_omega(x; operation=real), "area")
    #plot_n_save(x -> plot_soliton_z_lineshapes(x;nslices=4), "soliton_lineshapes")
    #plot_n_save(x -> plot_superimposed_lineshapes(x; constdim=(3,result.cfg.Ny÷2+1), nslices=10), "superimposed_slices")
    #plot_n_save(animate_field_2d, "anim")

    #plot_n_save(plot_2d_tz_ysum_heatmap, "testing")
    #plot_n_save(x -> plot_2d_z_ysum_line(x; operation=abs2), "testing")
    #plot_n_save(x -> plot_2d_z_ysum_line(x; operation=real), "testing")
    plot_n_save(plots_for_jevon, "testing")

    # Saving jld2 data
    #save_data(result, elapsed, "prop_2d")
end