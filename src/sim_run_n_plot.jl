
const PULSE_AREA_RTOL = 1e-8
const PULSE_AREA_ATOL = 1e-10

_area_is(area, target) = isapprox(area, target; rtol=PULSE_AREA_RTOL, atol=PULSE_AREA_ATOL)

"""
    classify_simulation(cfg) -> Symbol

Infer the physical test represented by `cfg`. The classification uses the pulse
sequence, coupling (`alpha`), diffraction (`beta`), and simulated spatial axes.
It throws an `ArgumentError` when the parameters do not describe one of the
supported tests.
"""
function classify_simulation(cfg)
    npulses = length(cfg.pulses)
    npulses in (1, 2) ||
        throw(ArgumentError("expected one or two pulses, got $npulses"))

    cfg.alpha >= 0 || throw(ArgumentError("alpha must be non-negative"))
    cfg.beta >= 0 || throw(ArgumentError("beta must be non-negative"))

    areas = [pulse[1].area for pulse in cfg.pulses]
    all(>(0), areas) || throw(ArgumentError("pulse areas must be positive; got $(areas)"))

    if npulses == 2
        is_echo = _area_is(areas[1], pi / 2) && _area_is(areas[2], pi)
        is_echo || throw(ArgumentError(
            "a two-pulse echo test requires pulse areas pi/2 and pi; got $(areas)",
        ))

        if cfg.Nz > 1 && iszero(cfg.alpha)
            throw(ArgumentError(
                "a propagating echo test requires alpha > 0 so the echo field can be emitted",
            ))
        end

        return :echo
    end

    # With no propagation axis, only the local atomic response is being tested.
    cfg.Nz == 1 && return :zero_dimensional

    has_medium_interaction = !iszero(cfg.alpha)
    has_diffraction = cfg.Ny > 1 && !iszero(cfg.beta)

    if !has_medium_interaction
        has_diffraction || throw(ArgumentError(
            "the configuration has neither medium interaction nor transverse diffraction",
        ))
        return :diffraction
    end

    medium_test = areas[1] > pi ? :sit : :attenuation

    if has_diffraction
        return medium_test == :sit ? :sit_with_diffraction : :attenuation_with_diffraction
    end

    return medium_test
end

function run_2d_simulation(cfg)
    simulation_test = classify_simulation(cfg)
    println("Detected simulation test: $(simulation_test)")

    # RUN
    @time result = run_2d_propagation(cfg)

    # Setting info for saving into file
    timestamp = Dates.format(now(), dateformat"yyyymmdd-HHMMSS-sss")
    parameter_info = "_test=$(simulation_test)_Nt=$(cfg.Nt)_Nd=$(cfg.Nd)_Nz=$(cfg.Nz)_Ny=$(cfg.Ny)_dwidth=$(cfg.d_width)_alpha=$(cfg.alpha)_beta=$(cfg.beta)_pulsecount=$(length(cfg.pulses))_"
    for pulse in cfg.pulses
        parameter_info *= "_area=$(pulse[1].area)_width=$(pulse[1].width)_"
    end


    # Setting plot directory and helper function
    plot_output_dir = joinpath(dirname(@__DIR__), "plots", "results")
    mkpath(plot_output_dir)
    plot_n_save(func, name) = save_plot(result, func, plot_output_dir, name; parameter_info=parameter_info, timestamp=timestamp)


    # Plotting
    if cfg.Ny > 1 && cfg.Nz > 1
        if simulation_test == :echo
            plot_n_save(plots_for_jevon_2sd_echo, "testing/finalstretch")
        elseif simulation_test == :diffraction
            plot_n_save(plot_fitted_diffraction, "testing/finalstretch")
        elseif simulation_test in (
            :sit,
            :attenuation,
            :sit_with_diffraction,
            :attenuation_with_diffraction,
        )
            plot_n_save(plots_for_jevon_2sd_SIT, "testing/finalstretch")

            simulation_test == :sit &&
                plot_n_save(plot_fitted_SIT, "testing/finalstretch")
            simulation_test == :attenuation &&
                plot_n_save(plot_fitted_attenuation, "testing/finalstretch")
        else
            error("No 2sd plots are defined for $(simulation_test)")
        end
    elseif cfg.Ny == 1 && cfg.Nz > 1 # 1sd
        if simulation_test == :echo
            plot_n_save(plots_for_jevon_1sd_echo, "testing/finalstretch")
        elseif simulation_test in (:sit, :attenuation)
            plot_n_save(plots_for_jevon_1sd_SIT, "testing/finalstretch")

            fit_plot = simulation_test == :sit ? plot_fitted_SIT : plot_fitted_attenuation
            plot_n_save(fit_plot, "testing/finalstretch")
        else
            error("No 1sd plots are defined for $(simulation_test)")
        end
    elseif cfg.Ny == 1 && cfg.Nz == 1 # 0sd
        plot_n_save(plots_for_jevon_0sd, "testing/finalstretch")
    else
        error("Unsupported spatial grid: Ny=$(cfg.Ny), Nz=$(cfg.Nz)")
    end

    # Saving jld2 data
    #save_data(result, elapsed, "prop_2d")

end
