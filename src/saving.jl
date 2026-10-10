using JLD2

function save_result(filename, result, elapsed_seconds)
    jldsave(
        filename;
        config=Dict(key => getfield(result.cfg, key) for key in fieldnames(typeof(result.cfg))),
        elapsed_seconds=elapsed_seconds,
        detunings=collect(result.detunings),
        time_vec=collect(result.time_vec),
        z_vec=collect(result.z_vec),
        Omega=result.Omega,
        P=result.P,
    )
end

function make_param_info(cfg)
    parameter_info = "_Nt=$(cfg.Nt)_Nd=$(cfg.Nd)_Nz=$(cfg.Nz)_Ny=$(cfg.Ny)_dwidth=$(cfg.d_width)_alpha=$(cfg.alpha)_beta=$(cfg.beta)_pulsecount=$(length(cfg.pulses))_"
    for pulse in cfg.pulses
        parameter_info *= "_area=$(pulse[1].area)_width=$(pulse[1].width)_"
    end
    return parameter_info
end

function save_plot(result, plot, plot_output_dir, subdir_name; parameter_info="placeholder", timestamp="placeholder")
    (parameter_info == "placeholder") && (parameter_info=make_param_info(result.cfg))
    (timestamp == "placeholder") && (timestamp=Dates.format(now(), dateformat"yyyymmdd-HHMMSS-sss"))

    fig = plot(result)
    (isnothing(fig)) && (return nothing)

    animbool = typeof(fig) <: Animation
    extension = (animbool) ? ".gif" : ".png"
    type = (animbool) ? "animation" : "plot"

    output_subdir = joinpath(plot_output_dir, subdir_name)
    mkpath(output_subdir)
    randomnumb = rand(1:1000)
    plot_path = joinpath(output_subdir, parameter_info * "_$(timestamp)_r$(randomnumb)" * extension)
    (animbool) ? gif(fig, plot_path) : savefig(fig, plot_path)

    path_elems = split(plot_path, "/")
    path_from_projroot = path_elems[end-2] * "/" * path_elems[end-1] * "/"
    println("Saved " * type * " to directory .../PROJECT_ROOT/" * path_from_projroot)
end

function save_data(result, elapsed, subdir_name; parameter_info="placeholder", timestamp="placeholder")
    (parameter_info == "placeholder") && (parameter_info=make_param_info(result.cfg);)
    (timestamp == "placeholder") && (timestamp=Dates.format(now(), dateformat"yyyymmdd-HHMMSS-sss");)

    data_output_dir = joinpath(dirname(@__DIR__), "data", subdir_name)
    mkpath(data_output_dir)
    data_path = joinpath(data_output_dir, parameter_info * "_$(timestamp)_" * ".jld2")
    save_result(data_path, result, elapsed)

    path_elems = split(data_path, "/")
    path_from_projroot = path_elems[end-2] * path_elems[end-1]

    println("Saved jld2 data to directory .../PROJECT_ROOT/" * path_from_projroot)
end;