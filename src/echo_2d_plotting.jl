using Plots, Measures


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

    output_subdir = joinpath(plot_output_dir, subdir_name);
    mkpath(output_subdir);
    randomnumb = rand(1:1000)
    plot_path = joinpath(output_subdir, parameter_info * "_$(timestamp)_r$(randomnumb)" * extension);
    (animbool) ? gif(fig, plot_path) : savefig(fig, plot_path);

    path_elems = split(plot_path, "/")
    path_from_projroot = path_elems[end-2] * "/" * path_elems[end-1] * "/"
    println("Saved " * type * " to directory .../PROJECT_ROOT/" * path_from_projroot);
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
    (result.cfg.Nz==1) && (return plot(result.time_vec, operation.(result.P[:, 1, 1]); label=false, xticks=false, yticks=false, c=:black, title="Polarisation"))
    nslices = clamp(nslices, 1, result.cfg.Nz)

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


function plot_superimposed_lineshapes(
    result;
    nslices=10,
    xaxisdim=1,
    constdim=(dim=3, value=1),
    operation=abs2,
    line_palette=:viridis,
)
    const_axis, const_index = constdim
    axis_vectors = (result.time_vec, result.z_vec, result.y_vec)
    axis_names = ("t", "z", "y")

    xaxisdim in 1:3 || throw(ArgumentError("xaxisdim must be 1, 2, or 3"))
    const_axis in 1:3 || throw(ArgumentError("constdim[1] must be 1, 2, or 3"))
    xaxisdim != const_axis ||
        throw(ArgumentError("xaxisdim and the constant dimension must be different"))

    # selectdim removes the constant axis while preserving the order of the
    # other two axes. For example, fixing y produces a (time, z) matrix.
    omega_slice = operation.(selectdim(result.Omega, const_axis, const_index))
    remaining_dims = filter(!=(const_axis), collect(1:3))
    matrix_xdim = only(findall(==(xaxisdim), remaining_dims))
    slice_dim = only(filter(!=(xaxisdim), remaining_dims))

    return plot_superimposed_slices(
        axis_vectors[xaxisdim],
        omega_slice;
        nslices,
        xdim=matrix_xdim,
        xlabel=axis_names[xaxisdim],
        ylabel="$(nameof(operation))(Ω)",
        slice_values=axis_vectors[slice_dim],
        slice_name=axis_names[slice_dim],
        line_palette,
        title="$(nameof(operation))(Ω) against $(axis_names[xaxisdim]) at constant $(axis_names[constdim[1]])=$(round(axis_vectors[constdim[1]][constdim[2]]; sigdigits=2))"
    )
end

function plot_superimposed_slices(
    xarray,
    yarray;
    nslices=10,
    xdim=1,
    xlabel="",
    ylabel="",
    slice_values=nothing,
    slice_name="slice",
    line_palette=:viridis,
    title=false
)
    ndims(yarray) == 2 || throw(ArgumentError("yarray must be two-dimensional"))
    xdim in 1:2 || throw(ArgumentError("xdim must be 1 or 2"))
    nslices >= 1 || throw(ArgumentError("nslices must be at least 1"))
    length(xarray) == size(yarray, xdim) || throw(
        DimensionMismatch(
            "x-axis has $(length(xarray)) points, but dimension $xdim of yarray has $(size(yarray, xdim))",
        ),
    )

    series_dim = 3 - xdim
    nseries = size(yarray, series_dim)
    slice_values !== nothing &&
        length(slice_values) != nseries &&
        throw(
            DimensionMismatch(
                "slice_values has $(length(slice_values)) points, but yarray contains $nseries slices",
            ),
        )

    # Include both ends and never request more lines than the data contains.
    slice_indices = unique(round.(Int, range(1, nseries; length=min(nslices, nseries))))
    line_colors = palette(line_palette, length(slice_indices))
    fig = plot(; xlabel, ylabel, legend=:topright, title)
    for (line_number, index) in enumerate(slice_indices)
        yseries = vec(selectdim(yarray, series_dim, index))
        is_first = line_number == firstindex(slice_indices)
        is_last = line_number == lastindex(slice_indices)
        label = if slice_values === nothing
            is_first ? "first $slice_name" : (is_last ? "last $slice_name" : false)
        else
            value = round(slice_values[index]; sigdigits=4)
            is_first ? "first: $slice_name = $value" :
            (is_last ? "last: $slice_name = $value" : false)
        end
        plot!(
            fig,
            xarray,
            yseries;
            color=line_colors[line_number],
            linewidth=2,
            label,
        )
    end

    return fig
end


function plot_sum_omega(result; operation=abs2)
    (result.cfg.Nz==1) && (return nothing)
    total = vec(sum(sum(operation.(result.Omega), dims=3), dims=1))
    title_str = "Sum of " * string(nameof(operation)) * " Omega"

    plot(result.z_vec, total, title=title_str, label=false, xlabel="z", ylims=(0.8*minimum(total), 1.2*maximum(total)))
end


#---------------------------

function plot_2d_tz_ysum_heatmap(result; operation=abs2, title=:default)
    Omega = operation.(dropdims(sum(result.Omega; dims=3); dims=3))
    Omega ./= maximum(Omega)
    t_vec, z_vec = result.time_vec, result.z_vec
    title = (title==:default) ? "sum $(nameof(operation)) Omega over y" : title
    return heatmap(t_vec, z_vec, transpose(Omega[:, :]); 
    c=:viridis, xlabel="t", ylabel="z", title=title, cbar=false)
end

function plot_2d_z_ysum_line(result; operation=abs2, title=:default, ylim_zoom_bool=false)
    Omega = operation.(dropdims(sum(result.Omega; dims=3); dims=3))
    z_vec = result.z_vec

    Omega_z = vec(sum(Omega; dims=1))
    Omega_z ./= Omega_z[begin]

    min, max = extrema(Omega_z)
    ymin = (all(>(0), Omega_z)) ? 0.0 : min-0.05abs(min)
    ymin = (ylim_zoom_bool) ? min-0.05abs(min) : ymin;
    ymax = max+0.05abs(max)
    ylims = (ymin, ymax)
    x_axis = z_vec
    xlims = extrema(x_axis)

    title = (title==:default) ? "sum $(nameof(operation)) Omega over y and t" : title

    return plot(x_axis, Omega_z; xlabel="z", title=title, legend=false, ylims=ylims, xlims=xlims)
end

function plot_2d_superimposed_z_ysum_line(result; operation=abs2, nslices=10, title=:default)
    Omega = operation.(dropdims(sum(result.Omega; dims=3); dims=3))
    Omega ./= maximum(Omega)
    
    t_vec, z_vec = result.time_vec, result.z_vec
    title = (title==:default) ? "sum $(nameof(operation)) Omega over y" : title
    min, max = extrema(Omega)
    ylims = extrema(Omega)
    x_axis = t_vec
    xlims = extrema(x_axis) #watch it

    colgrad = cgrad([:blue, :yellow, :red])
    colors = colgrad[range(0, 1, length=nslices)]
    fig = plot(; xlabel="t", title=title, ylims=ylims, xlims=xlims)

    firstlast = i -> (i == 1) ? "first" : ((i==nslices) ? "last" : ((i==nslices÷2+1) ? "middle" : false))

    delta = floor(Int, length(z_vec) / nslices)


    for i in 1:nslices
        plot!(x_axis, Omega[:, i*delta-(delta-1)], c=colors[i], label=firstlast(i))
    end

    return fig
end


function plot_2d_z_yslices_line(result; operation=abs2, nslices=10, title=:default)
    Omega = operation.(result.Omega[:, :, (end÷2+1):end])
    z_vec = result.z_vec
    title = (title==:default) ? "sum $(nameof(operation)) Omega and t" : title
    Omega_z = dropdims(sum(Omega; dims=1); dims=1)
    Omega_z ./= maximum(Omega_z)
    min, max = extrema(Omega_z)
    ylims = (0.0, max+0.05abs(max))
    x_axis = z_vec
    xlims = extrema(x_axis)

    colgrad = cgrad([:blue, :yellow, :red])
    colors = colgrad[range(0, 1, length=nslices)]
    fig = plot(; xlabel="z", title=title, ylims=ylims, xlims=xlims, legend=false)

    delta = floor(Int, length(result.y_vec) / 2nslices)


    for i in 1:nslices
        plot!(x_axis, Omega_z[:, i*delta-(delta-1)], c=colors[i])
    end

    return fig
end

function plot_2d_zy_tslices_heatmap(result; operation=abs2, nslices=3, title=:default)
    Omega = operation.(result.Omega)
    Omega ./= maximum(Omega)
    y_vec, z_vec = result.y_vec, result.z_vec
    
    fig_vec = Array{Plots.Plot}(undef,nslices)
    index = i -> ceil(Int, (length(result.time_vec) / nslices) * i)

    for i in 1:nslices
        fig_vec[i] = heatmap(y_vec, z_vec, transpose(Omega[index(i),:,:]); c=:viridis, xlabel="y", ylabel="z")
    end

    return plot(fig_vec...; layout=(nslices, 1), size=(600,200*nslices))
end

function plot_2d_zy_pulseprofile_heatmap(result; operation=abs2, title=:default)
    Omega = operation.(result.Omega)
    Omega ./= maximum(Omega)
    y_vec, z_vec = result.y_vec, result.z_vec

    tindex = argmax(Omega[:,end÷2+1,end÷2+1])
    Omega_tslice = Omega[tindex,:,:]
    Omega_tslice ./= maximum(Omega_tslice)

    title = (title==:default) ? "$(nameof(operation)) Omega time slice" : title

    return heatmap(y_vec, z_vec, Omega_tslice; c=:viridis, xlabel="y", ylabel="z", title=title, cbar=false)
end

function plot_2sd_spatialslices_superimposed(result; nslices=4)
    Omega = abs2.(result.Omega)
    Omega ./= maximum(Omega)

    halfOmegas = filter(x -> (x<=0.501 && x>=0.499), Omega)

    
    for i in 1:nslices
        plot(;seriestype=:scatter)
    end
end

function plots_for_jevon_2sd(result)
    heatmap = plot_2d_tz_ysum_heatmap(result; operation=abs2, title="Pulse intensity")
    #area = plot_2d_z_ysum_line(result; operation=real, title="Total pulse area")
    pulse_profile_heatmap = plot_2d_zy_pulseprofile_heatmap(result; operation=abs2, title="Pulse intensity profile time slice")
    energy = plot_2d_z_ysum_line(result; operation=abs2, title="Total pulse energy")
    superimposed = plot_2d_superimposed_z_ysum_line(result; operation=abs2, nslices=20, title="Pulse intensity slices in z")
    superimposed_area_yslices = plot_2d_z_yslices_line(result; operation=real, nslices=30, title="Pulse area slices in y against z")
    superimposed_intensity_yslices = plot_2d_z_yslices_line(result; operation=abs2, nslices=30, title="Pulse intensity slices in y against z")

    return plot(
        heatmap, pulse_profile_heatmap, superimposed_area_yslices, 
        superimposed, energy, superimposed_intensity_yslices, 
        layout=(2, 3), size=(1600, 1000), margins=3mm, framestyle=:box, left_margin=5mm, bottom_margin=5mm)
end


#-------------------------------------------

function plot_1d_tz_heatmap(result; operation=abs2, title=:default)
    Omega = operation.(dropdims(result.Omega; dims=3))
    Omega ./= maximum(Omega)
    t_vec, z_vec = result.time_vec, result.z_vec
    title = (title==:default) ? "$(nameof(operation)) Omega" : title
    return heatmap(t_vec, z_vec, transpose(Omega[:, :]); 
    c=:viridis, xlabel="t", ylabel="z", title=title)
end

function plot_1d_superimposed_lines_against_t(result; operation=abs2, nslices=10, title=:default, zoom=false)
    Omega = operation.(dropdims(result.Omega;dims=3))
    Omega ./= maximum(Omega)
    
    t_vec, z_vec = result.time_vec, result.z_vec
    title = (title==:default) ? "$(nameof(operation)) Omega" : title
    min, max = extrema(Omega)
    ylims = extrema(Omega)
    x_axis = t_vec
    xlims = extrema(x_axis) #watch it

    if zoom
        startindex = findfirst(x -> x>0.01, Omega[:,begin])
        endindex = length(t_vec) - findfirst(x -> x>0.01, reverse(Omega[:,end]))

        x_axis = x_axis[startindex:endindex]
        Omega = Omega[startindex:endindex,:]

        xlims = extrema(x_axis)
    end

    colgrad = cgrad([:blue, :yellow, :red])
    colors = colgrad[range(0, 1, length=nslices)]
    fig = plot(; xlabel="t", title=title, ylims=ylims, xlims=xlims)

    firstlast = i -> (i == 1) ? "first" : ((i==nslices) ? "last" : ((i==nslices÷2+1) ? "middle" : false))

    delta = floor(Int, length(z_vec) / nslices)

    for i in 1:nslices
        plot!(x_axis, Omega[:, i*delta-(delta-1)], c=colors[i], label=firstlast(i))
    end

    return fig
end

function plot_1d_superimposed_lines_against_z(result; operation=abs2, nslices=10, title=:default)
    Omega = operation.(dropdims(result.Omega; dims=3))
    Omega ./= maximum(Omega)
    
    t_vec, z_vec = result.time_vec, result.z_vec
    title = (title==:default) ? "$(nameof(operation)) Omega" : title
    min, max = extrema(Omega)
    ylims = (0.0, max+0.05abs(max))
    x_axis = z_vec
    xlims = extrema(x_axis) #watch it

    colgrad = cgrad([:blue, :yellow, :red])
    colors = colgrad[range(0, 1, length=nslices)]
    fig = plot(; xlabel="z", title=title, ylims=ylims, xlims=xlims)

    firstlast = i -> (i == 1) ? "first" : ((i==nslices) ? "last" : ((i==nslices÷2+1) ? "middle" : false))

    delta = floor(Int, length(t_vec) / nslices)


    for i in 1:nslices
        plot!(x_axis, Omega[i*delta-(delta-1), :], c=colors[i], label=firstlast(i))
    end

    return fig
end


function plots_for_jevon_1sd(result)
    heatmap = plot_1d_tz_heatmap(result; operation=abs2, title="Pulse intensity")
    superimposed_zlines_against_t = plot_1d_superimposed_lines_against_t(result; title="Pulse intensity slices in \$z\$", zoom=true, nslices=10)
    energy = plot_2d_z_ysum_line(result; operation=abs2, title="Sum of pulse intensity across \$t\$", ylim_zoom_bool=true)
    area = plot_2d_z_ysum_line(result; operation=real, title="Pulse area", ylim_zoom_bool=true)
    #superimposed_tlines_against_z = plot_1d_superimposed_lines_against_z(result; operation=abs2, nslices=10, title="Pulse intensity slices in \$t\$")

    return plot(
        heatmap, energy, 
        superimposed_zlines_against_t, area,
        layout=(2, 2), size=(800, 500), margins=2mm, framestyle=:box, left_margin=4mm, bottom_margin=4mm)
end

#------------------------

function plot_0sd_polarisation(result)
    P = abs2.(result.P)

    t_vec = result.time_vec
    #title = (title==:default) ? "" : title
    fig = plot(; ylabel="\$|\\mathcal{P}|^2\$", title="Squared magnitude of the polarisation density against time")

    plot!(t_vec, P; label=false)

    return fig
end

function plot_0sd_pol_complex_parts(result)
    Pu = 2*real(result.P)
    Pv = -2*imag(result.P)

    t_vec = result.time_vec

    fig = plot(; title="Ensemble averaged Bloch components against time")

    plot!(t_vec, Pu; label="averaged u")
    plot!(t_vec, Pv; label="averaged v")

    return fig
end

function plots_for_jevon_0sd(result)
    pol_abs2 = plot_0sd_polarisation(result)
    pol_real_n_imag = plot_0sd_pol_complex_parts(result)


    return plot(
        pol_abs2, pol_real_n_imag; layout=(2,1), size=(800,600)
    )
end