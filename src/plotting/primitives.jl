"""
    _unit_scale_ylims(values; padding_fraction=0.05)

Return plot limits spanning both the data and `[0, 1]`, with fractional
padding above and, when values are negative, below the span.
"""
function _unit_scale_ylims(values; padding_fraction=0.05)
    lower = min(0.0, minimum(values))
    upper = max(1.0, maximum(values))
    span = upper - lower
    padded_lower = lower < 0 ? lower - padding_fraction * span : 0.0
    padded_upper = upper + padding_fraction * span
    return (padded_lower, padded_upper)
end

"""
    _plot_tz_heatmap(time, z, field; title)

Return a heatmap of `field[time, z]` with time on the horizontal axis
and longitudinal position on the vertical axis, using the shared style.
"""
function _plot_tz_heatmap(time, z, field; title)
    return heatmap(time, z, transpose(field);
        c=:viridis, xlabel="\$t\$", ylabel="\$z\$", title=title,
        cbar=false, size=PLOT_PANEL_SIZE)
end

"""
    _plot_xlims(x, field, zoom)

Return the extrema of coordinate vector `x`, or, when `zoom` is true,
the coordinates bounding above-threshold rows of the field matrix
according to [`zoomingin`](@ref).
"""
function _plot_xlims(x, field, zoom)
    zoom || return extrema(x)
    startindex, endindex = zoomingin(field)
    return (x[startindex], x[endindex])
end

"""
    _plot_temporal_slices(field, time, indices; title, zoom, ylabel, size=PLOT_PANEL_SIZE)

Return overlaid curves of `field[:, index]` against `time` for the
selected column `indices`. Use shared slice colours and legend labels;
`zoom=true` restricts time limits to the field above the default threshold.
"""
function _plot_temporal_slices(field, time, indices; title, zoom, ylabel, size=PLOT_PANEL_SIZE)
    xlims = _plot_xlims(time, field, zoom)
    colors = _slice_colors(length(indices))
    fig = plot(; xlabel="\$t\$", ylabel=ylabel, title=title, xlims=xlims, size=size)

    for (i, field_index) in enumerate(indices)
        plot!(time, field[:, field_index]; color=colors[i], label=_slice_label(i, length(indices)))
    end

    return fig
end
