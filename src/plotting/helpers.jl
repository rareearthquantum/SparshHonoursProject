using Plots, Measures

const PLOT_PANEL_SIZE = (800, 600)
const PLOT_LEFT_MARGIN = 2mm
const PLOT_BOTTOM_MARGIN = 2mm

"""
    _plot_margins(layout)

Return `(left_margin, bottom_margin)` for the requested subplot layout,
falling back to the shared defaults for other layouts.
"""
function _plot_margins(layout)
    layout == (1, 2) && return (5mm, 3mm)
    layout == (1, 3) && return (3mm, 4mm)
    layout == (2, 2) && return (2mm, 1mm)
    return (PLOT_LEFT_MARGIN, PLOT_BOTTOM_MARGIN)
end

"""
    _normalise_by_max!(values)

Divide every entry of `values` by its maximum in place and return the
same array. Requires a nonzero maximum for finite normalised values.
"""
function _normalise_by_max!(values)
    values ./= maximum(values)
    return values
end

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
    _field_2d_summed_over_y(result, operation)

Apply `operation` elementwise to `result.Omega[time, z, y]` and sum
over the transverse axis. Return a time-by-z matrix; no grid-spacing
factor or normalisation is applied.
"""
_field_2d_summed_over_y(result, operation) =
    dropdims(sum(operation.(result.Omega); dims=3); dims=3)

"""
    _field_1d(result, operation)

Remove the singleton transverse axis of `result.Omega` and apply
`operation` elementwise, returning a time-by-z matrix. Requires `Ny == 1`.
"""
_field_1d(result, operation) = operation.(dropdims(result.Omega; dims=3))

"""
    _slice_colors(count)

Return `count` evenly spaced colours from the blue–yellow–red gradient.
"""
_slice_colors(count) = cgrad([:blue, :yellow, :red])[range(0, 1, length=count)]

"""
    _slice_label(index, count)

Return a legend label for the first, last, or middle slice; return
`false` to omit other slices from the legend. `index` is the slice ordinal.
"""
function _slice_label(index, count)
    index == 1 && return "first"
    index == count && return "last"
    index == count ÷ 2 + 1 && return "middle"
    return false
end

"""
    _add_panel_titles!(plots)

Replace subplot titles with letter labels `(a)`, `(b)`, and so on,
positioned on the left. Mutate each panel and return `plots`.
"""
function _add_panel_titles!(plots)
    for (index, panel) in enumerate(plots)
        panel_label = "($(Char(Int('a') + index - 1)))"
        title!(panel, panel_label; titlelocation=:left)
    end
    return plots
end

"""
    _compose_plots(plots...; layout)

Label the supplied panels and return a combined plot with the requested
`(rows, columns)` layout, shared margins, and scaled figure height.
The input panel titles are mutated.
"""
function _compose_plots(plots...; layout)
    _add_panel_titles!(plots)
    rows, columns = layout
    height = round(Int, PLOT_PANEL_SIZE[2] * rows / columns)
    size = (PLOT_PANEL_SIZE[1], height)
    left_margin, bottom_margin = _plot_margins(layout)
    return plot(plots...; layout=layout, size=size,
        left_margin=left_margin, bottom_margin=bottom_margin)
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

"""
    zoomingin(vector; threshold=0.001)

Return inclusive row bounds around entries at least `threshold` in any
column of the matrix `vector`, padded by one row within array bounds.
Requires at least one above-threshold entry for meaningful ordered bounds.
"""
function zoomingin(vector; threshold=0.001)

    first_above = lastindex(vector[:, 1])
    last_above = firstindex(vector[:, 1])

    for i in 1:length(vector[1, :])
        viewvec = @view vector[:, i]
        temp = findfirst(>=(threshold), viewvec)
        if temp !== nothing && temp < first_above
            first_above = temp
        end

        temp = findlast(>=(threshold), viewvec)
        if temp !== nothing && temp > last_above
            last_above = temp
        end
    end

    startindex = max(first_above - 1, firstindex(vector[:, 1]))
    endindex = min(last_above + 1, lastindex(vector[:, 1]))

    return (startindex, endindex)
end
