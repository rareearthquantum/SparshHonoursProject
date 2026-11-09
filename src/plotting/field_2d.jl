"""
    plot_2d_tz_ysum_heatmap(result; operation=abs2, title=:default)

Return a time–z heatmap of `operation.(Omega)` summed over transverse
positions and normalised by its maximum. `title=:default` derives a title.
"""
function plot_2d_tz_ysum_heatmap(result; operation=abs2, title=:default)
    field_values = _field_2d_summed_over_y(result, operation)
    _normalise_by_max!(field_values)
    title = (title==:default) ? "sum $(nameof(operation)) Omega over y" : title
    return _plot_tz_heatmap(result.t_vec, result.z_vec, field_values; title=title)
end

"""
    plot_2d_z_ysum_line(result; operation=abs2, title=:default, size=PLOT_PANEL_SIZE, ylabel="", label=false, unit_scale=false)

Return the field summed over time and transverse positions versus z,
after applying `operation` and dividing by the input-plane sum.
`unit_scale=true` extends the vertical limits to include `[0, 1]`.
The input sum must be nonzero; sums omit grid-spacing factors.
"""
function plot_2d_z_ysum_line(result; operation=abs2, title=:default,
        size=PLOT_PANEL_SIZE, ylabel="", label=false, unit_scale=false)
    field_values = _field_2d_summed_over_y(result, operation)
    z_vec = result.z_vec

    longitudinal_field_values = vec(sum(field_values; dims=1))
    longitudinal_field_values ./= longitudinal_field_values[begin]

    x_axis = z_vec
    xlims = extrema(x_axis)

    title = (title==:default) ? "sum $(nameof(operation)) Omega over y and t" : title

    fig = plot(x_axis, longitudinal_field_values; xlabel="\$z\$", ylabel=ylabel,
        title=title, label=label, xlims=xlims, size=size)
    unit_scale && ylims!(fig, _unit_scale_ylims(longitudinal_field_values))
    return fig
end

"""
    plot_2d_superimposed_z_ysum_line(result; operation=abs2, nslices=10, title=:default, zoom=true, ylabel="")

Return temporal profiles of the transversely summed transformed field
at `nslices` longitudinal indices, using global maximum normalisation.
`zoom=true` crops time to the above-threshold region. Choose
`1 <= nslices <= length(result.z_vec)` for distinct samples.
"""
function plot_2d_superimposed_z_ysum_line(result; operation=abs2, nslices=10, title=:default, zoom=true, ylabel="")
    field_values = _field_2d_summed_over_y(result, operation)
    _normalise_by_max!(field_values)

    title = (title==:default) ? "sum $(nameof(operation)) Omega over y" : title
    delta = floor(Int, length(result.z_vec) / nslices)
    indices = 1 .+ (0:nslices-1) .* delta
    return _plot_temporal_slices(field_values, result.t_vec, indices;
        title=title, zoom=zoom, ylabel=ylabel)
end

"""
    plot_2d_z_yslices_line(result; operation=abs2, nslices=10, title=:default, ylabel="")

Return time-summed longitudinal curves at `nslices` transverse indices
from the second half of the transverse grid, after applying `operation`
and shared maximum normalisation. Choose `nslices` no larger than half
the transverse sample count for distinct slices.
"""
function plot_2d_z_yslices_line(result; operation=abs2, nslices=10, title=:default, ylabel="")
    field_values = operation.(result.Omega[:, :, (end÷2+1):end])
    z_vec = result.z_vec
    title = (title==:default) ? "sum $(nameof(operation)) Omega and t" : title
    longitudinal_field_values = dropdims(sum(field_values; dims=1); dims=1)
    _normalise_by_max!(longitudinal_field_values)
    x_axis = z_vec
    xlims = extrema(x_axis)

    colors = _slice_colors(nslices)
    fig = plot(; xlabel="\$z\$", ylabel=ylabel, title=title, xlims=xlims,
        legend=false, size=PLOT_PANEL_SIZE)

    delta = floor(Int, length(result.y_vec) / 2nslices)


    for i in 1:nslices
        plot!(x_axis, longitudinal_field_values[:, i*delta-(delta-1)], c=colors[i])
    end

    return fig
end

"""
    plot_2d_zy_pulseprofile_heatmap(result; operation=abs2, title=:default)

Return a transverse–longitudinal heatmap at the time of peak transformed
field at the central z and y indices. Apply `operation` and normalise
the selected time slice by its maximum.
"""
function plot_2d_zy_pulseprofile_heatmap(result; operation=abs2, title=:default)
    field_values = operation.(result.Omega)
    _normalise_by_max!(field_values)
    y_vec, z_vec = result.y_vec, result.z_vec

    tindex = argmax(field_values[:, end÷2+1, end÷2+1])
    field_time_slice = field_values[tindex, :, :]
    _normalise_by_max!(field_time_slice)

    title = (title==:default) ? "$(nameof(operation)) Omega time slice" : title

    return heatmap(y_vec, z_vec, field_time_slice; c=:viridis,
        xlabel="\$y\$", ylabel="\$z\$", title=title,
        cbar=false, size=PLOT_PANEL_SIZE)
end

"""
    plot_2d_superimposed_lines_against_t(result; operation=abs2, nslices=10, title="", zoom=false, ylabel="")

Return overlaid temporal profiles of `operation.(Omega)` summed over y
at `nslices` longitudinal indices, normalised by their global maximum.
`zoom=true` crops the time range. Choose `nslices` no larger than the
longitudinal sample count for distinct slices.
"""
function plot_2d_superimposed_lines_against_t(result; operation=abs2, nslices=10, title="", zoom=false, ylabel="")
    field_values = _field_2d_summed_over_y(result, operation)
    _normalise_by_max!(field_values)

    title = (title==:default) ? "$(nameof(operation)) Omega" : title
    delta = floor(Int, length(result.z_vec) / nslices)
    indices = 1 .+ (0:nslices-1) .* delta
    return _plot_temporal_slices(field_values, result.t_vec, indices;
        title=title, zoom=zoom, ylabel=ylabel)
end

"""
    plot_2d_superimposed_lines_against_y(result; operation=abs2, nslices=10, title=:default, size=PLOT_PANEL_SIZE, ylabel="")

Return transverse profiles at `nslices` evenly spaced longitudinal
indices at the middle time sample. Apply `operation` and normalise
by the maximum over that time slice; `size` and `ylabel` style the plot.
"""
function plot_2d_superimposed_lines_against_y(result; operation=abs2, nslices=10, title=:default, size=PLOT_PANEL_SIZE, ylabel="")
    field_values = operation.(result.Omega[end÷2+1, :, :])
    _normalise_by_max!(field_values)

    y_vec = result.y_vec
    z_vec = result.z_vec
    title = (title==:default) ? "$(nameof(operation)) Omega" : title
    x_axis = y_vec
    xlims = extrema(x_axis) #watch it

    colors = _slice_colors(nslices)
    fig = plot(; xlabel="\$y\$", ylabel=ylabel, title=title, xlims=xlims, size=size)

    indices = round.(Int, range(1, length(z_vec), length=nslices))

    for (i, j) in enumerate(indices)
        plot!(x_axis, field_values[j, :], c=colors[i], label=_slice_label(i, nslices))
    end

    return fig
end
