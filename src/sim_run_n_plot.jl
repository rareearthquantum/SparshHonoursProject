
function run_2d_simulation(cfg)
    # RUN
    @time result = run_2d_propagation(cfg)

    # Setting info for saving into file
    timestamp = Dates.format(now(), dateformat"yyyymmdd-HHMMSS-sss")
    parameter_info = "_Nt=$(cfg.Nt)_Nd=$(cfg.Nd)_Nz=$(cfg.Nz)_Ny=$(cfg.Ny)_dwidth=$(cfg.d_width)_alpha=$(cfg.alpha)_beta=$(cfg.beta)_pulsecount=$(length(cfg.pulses))_"
    for pulse in cfg.pulses
        parameter_info *= "_area=$(pulse[1].area)_width=$(pulse[1].width)_"
    end


    # Setting plot directory and helper function
    plot_output_dir = joinpath(dirname(@__DIR__), "plots", "results")
    mkpath(plot_output_dir)
    plot_n_save(func, name) = save_plot(result, func, plot_output_dir, name; parameter_info=parameter_info, timestamp=timestamp)


    #=
    # Plotting
    if (cfg.Ny > 1 && cfg.Nz > 1)
        
        if (length(cfg.pulses) == 1)
            plot_n_save(plots_for_jevon_2sd_SIT, "testing/2sd")
        elseif (length(cfg.pulses)==2)
            plot_n_save(plots_for_jevon_2sd_echo, "testing/2sd")
        else
            error("How many pulses?")
        end

    elseif (cfg.Ny == 1 && cfg.Nz > 1) #1sd

        if (length(cfg.pulses) == 1)
            plot_n_save(plots_for_jevon_1sd_SIT, "testing/1sd")
        elseif (length(cfg.pulses)==2)
            plot_n_save(plots_for_jevon_1sd_echo, "testing/1sd")
        else
            error("How many pulses?")
        end
        
    elseif (cfg.Ny==1 && cfg.Nz == 1) #0sd
        #clearly echo test
        plot_n_save(plots_for_jevon_0sd, "testing/0sd")
    else
        error("What kind of Ny and Nz do you have???")
    end

    # Saving jld2 data
    #save_data(result, elapsed, "prop_2d")
    =#

    plot_n_save(plot_fitted_diffraction, "testing/fits")

end