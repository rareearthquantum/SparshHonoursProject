function plot_1d_tz_heatmap(result; operation=abs2, title=:default)
    field_values = _field_1d(result, operation)
    _normalise_by_max!(field_values)
    title = (title==:default) ? "$(nameof(operation)) Omega" : title
    return _plot_tz_heatmap(result.t_vec, result.z_vec, field_values; title=title)
end

function plot_1d_superimposed_lines_against_t(result; operation=abs2, nslices=10, title=:default, zoom=false, size=PLOT_PANEL_SIZE, ylabel="")
    field_values = _field_1d(result, operation)
    _normalise_by_max!(field_values)

    title = (title==:default) ? "$(nameof(operation)) Omega" : title
    indices = round.(Int, range(1, length(result.z_vec), length=nslices))
    return _plot_temporal_slices(field_values, result.t_vec, indices;
        title=title, zoom=zoom, ylabel=ylabel, size=size)
end

function plot_1d_superimposed_lines_against_z(result; operation=abs2, nslices=10, title=:default)
    field_values = _field_1d(result, operation)
    _normalise_by_max!(field_values)

    t_vec, z_vec = result.t_vec, result.z_vec
    title = (title==:default) ? "$(nameof(operation)) Omega" : title
    x_axis = z_vec
    xlims = extrema(x_axis) #watch it

    colors = _slice_colors(nslices)
    fig = plot(; xlabel="\$z\$", title=title, xlims=xlims, size=PLOT_PANEL_SIZE)

    delta = floor(Int, length(t_vec) / nslices)


    for i in 1:nslices
        plot!(x_axis, field_values[i*delta-(delta-1), :], c=colors[i], label=_slice_label(i, nslices))
    end

    return fig
end
