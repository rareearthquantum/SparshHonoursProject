using Plots


function make_param_info(cfg)
    parameter_info = "_Nt=$(cfg.Nt)_Nd=$(cfg.Nd)_Nz=$(cfg.Nz)_Ny=$(cfg.Ny)_dwidth=$(cfg.d_width)_alpha=$(cfg.alpha)_beta=$(cfg.beta)_pulsecount=$(length(cfg.pulses))_"
    for pulse in cfg.pulses
        parameter_info *= "_area=$(pulse[1].area)_width=$(pulse[1].width)_"
    end
    return parameter_info
end

function save_plot(result, plot, plot_output_dir, subdir_name; parameter_info="placeholder", timestamp="placeholder")
    (parameter_info == "placeholder") && (parameter_info=make_param_info(result.cfg);)
    (timestamp == "placeholder") && (timestamp=Dates.format(now(), dateformat"yyyymmdd-HHMMSS-sss");)

    fig = plot(result)

    animbool = typeof(fig) <: Animation
    extension = (animbool) ? ".gif" : ".png"
    type = (animbool) ? "animation" : "plot"

    output_subdir = joinpath(plot_output_dir, subdir_name);
    mkpath(output_subdir);
    plot_path = joinpath(output_subdir, parameter_info * "_$(timestamp)_" * extension);
    (animbool) ? gif(fig, plot_path) : savefig(fig, plot_path);

    path_elems = split(plot_path, "/")
    path_from_projroot = path_elems[end-2] * "/" * path_elems[end-1] * "/"

    println("Saved " * type * " to directory .../PROJECT_ROOT/" * path_from_projroot);
end;

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




function animate_field_2d(result; filename="echo_2d.gif", fps=30, operation=abs2)
    field_intensity = operation.(result.Omega)
    z_vec = result.z_vec
    y_vec = result.y_vec
    t_vec = result.time_vec

    anim = @animate for i in eachindex(t_vec)
        heatmap(y_vec, z_vec, field_intensity[i, :, :];
            clabel="|Ω|²",
            ylabel="z",
            xlabel="y",
            clims=extrema(field_intensity),
            legend=false,
            title="Rabi frequency intensity at t = $(round(t_vec[i]; sigdigits=4))",
            c=:viridis
        )
    end

    return anim;
end


function plot_soliton_z_lineshapes(result; nslices=5, operation=abs2)
    nslices = clamp(nslices, 1, 10)

    Omega_tz = operation.(result.Omega[:, :, end÷2+1])
    t_vec = result.time_vec

    fig_vec = Array{Plots.Plot}(undef, nslices)
    fig_begin = plot(t_vec, Omega_tz[:, begin]; label=false, xticks=false, yticks=false, c=:black, ylims=extrema(Omega_tz))
    delta = cfg.Nz / nslices
    for i in 1:nslices
        zindex = floor(Int, delta*i)
        fig_vec[i] = plot(t_vec, Omega_tz[:, zindex]; label=false, xticks=false, yticks=false, c=:black, ylims=extrema(Omega_tz))
    end

    title_str = string(nameof(operation)) * " of Omega"

    fig = plot(fig_begin, fig_vec...; layout=(nslices+1, 1), size=(400, nslices*100))

    return fig
end

function plot_soliton_t_lineshapes(result; nslices=10, operation=abs2)
    nslices = clamp(nslices, 1, 10)

    Omega_tz = operation.(result.Omega[:, :, end÷2+1])
    z_vec = result.z_vec

    fig_vec = Array{Plots.Plot}(undef, nslices)
    fig_begin = plot(z_vec, Omega_tz[begin, :]; label=false, xticks=false, yticks=false, c=:black, ylims=extrema(Omega_tz))
    delta = cfg.Nt / nslices
    for i in 1:nslices
        tindex = floor(Int, delta*i)
        fig_vec[i] = plot(z_vec, Omega_tz[tindex, :]; label=false, xticks=false, yticks=false, c=:black, ylims=extrema(Omega_tz))
    end

    title_str = string(nameof(operation)) * " of Omega"

    fig = plot(fig_begin, fig_vec...; layout=(nslices+1, 1), size=(400, nslices*100))

    return fig
end



function plot_sum_omega(result; operation=abs2)
    total = vec(sum(sum(operation.(result.Omega), dims=3), dims=1))
    title_str = "Sum of " * string(nameof(operation)) * " Omega"

    plot(result.z_vec, total, title=title_str, label=false, xlabel="z", ylims=(0.8*minimum(total), 1.2*maximum(total)))
end


