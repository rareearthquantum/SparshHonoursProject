# Simulation setup and propagation loop.

function make_rotate_grid(detunings, time_vec)
    return [exp(im * δ * t) for δ in detunings, t in time_vec]
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

function atom_stage_params(Omega_col, rotate_vec, time_vec)
    # Currently all RK4 stages use the same field column Ω[:, j], matching
    # the uploaded working script. If you later add true half-step fields,
    # this is the place to pass different p1, p2, p3, p4 values.
    p = (Omega_col, rotate_vec, time_vec)
    return (p, p, p, p)
end

function compute_polarisation_column!(P_col, sigma_temp, Omega_col, rotate_grid, time_vec)
    Nd = size(rotate_grid, 1)
    fill!(P_col, 0)

    @inbounds for i in 1:Nd
        reset_atom_history!(sigma_temp)

        rotate_vec = @view rotate_grid[i, :]
        ps = atom_stage_params(Omega_col, rotate_vec, time_vec)

        rk4_staged!(atom_single!, sigma_temp, time_vec, ps)

        @views @. P_col += sigma_temp[1, :] * conj(rotate_vec)
    end

    @. P_col /= Nd

    return nothing
end

function run_propagation(cfg::EchoConfig=EchoConfig();
        omega_input=make_default_omega_input(cfg),
        compute_final_polarisation::Bool=true)

    detunings = make_detunings(cfg)
    time_vec = make_time_grid(cfg)
    z_vec = make_z_grid(cfg)
    dz = step(z_vec)

    Omega = initialise_field(cfg, time_vec, z_vec, omega_input)
    P = initialise_polarisation(time_vec, z_vec)
    sigma_temp = zeros(ComplexF64, 2, length(time_vec))

    rotate_grid = make_rotate_grid(detunings, time_vec)
    field_cache = @views AB2Cache(Omega[:, 1])

    @inbounds for j in 1:(length(z_vec)-1)
        @views compute_polarisation_column!(P[:, j], sigma_temp, Omega[:, j], rotate_grid, time_vec)

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

    # The original script left P[:, end] at zero because the final z column
    # is never needed to step Ω forward. Computing it avoids a misleading
    # zero stripe in the final polarisation heatmap.
    if compute_final_polarisation
        @views compute_polarisation_column!(P[:, end], sigma_temp, Omega[:, end], rotate_grid, time_vec)
    end

    return (
        cfg=cfg,
        detunings=detunings,
        time_vec=time_vec,
        z_vec=z_vec,
        Omega=Omega,
        P=P,
        rotate_grid=rotate_grid,
    )
end
