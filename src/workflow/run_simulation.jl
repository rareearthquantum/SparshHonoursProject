using Dates

"""
    run_2d_simulation(cfg)

Classify `cfg`, run the propagation solver with timing and progress
output, and save the selected research figures beneath `plots/results/`.
Filenames include the classification, parameters, and timestamp. Numerical
data saving is currently disabled. Returns `nothing` on success;
unsupported experiment configurations or plot selections raise an error.
"""
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
    plot_output_dir = joinpath(dirname(dirname(@__DIR__)), "plots", "results")
    mkpath(plot_output_dir)


    save_simulation_plots(result, simulation_test, plot_output_dir; parameter_info, timestamp)

    # Saving jld2 data
    #save_data(result, elapsed, "prop_2d")

end
