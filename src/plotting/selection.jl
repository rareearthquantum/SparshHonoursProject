"""Select and save the existing research plots for an already computed result."""
function save_simulation_plots(result, simulation_test, plot_output_dir; parameter_info, timestamp)
    cfg = result.cfg
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

            simulation_test in (:sit, :sit_with_diffraction) &&
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
    return nothing
end
