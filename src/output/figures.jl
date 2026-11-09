using Dates

function save_plot(result, plot, plot_output_dir, subdir_name; parameter_info="placeholder", timestamp="placeholder")
    (parameter_info == "placeholder") && (parameter_info=make_param_info(result.cfg))
    (timestamp == "placeholder") && (timestamp=Dates.format(now(), dateformat"yyyymmdd-HHMMSS-sss"))

    fig = plot(result)
    (isnothing(fig)) && (return nothing)

    output_subdir = joinpath(plot_output_dir, subdir_name)
    mkpath(output_subdir)
    randomnumb = rand(1:1000)
    plot_path = joinpath(output_subdir, parameter_info * "_$(timestamp)_r$(randomnumb).pdf")
    savefig(fig, plot_path)

    path_elems = split(plot_path, "/")
    path_from_projroot = path_elems[end-2] * "/" * path_elems[end-1] * "/"
    println("Saved plot to directory .../PROJECT_ROOT/" * path_from_projroot)
end
