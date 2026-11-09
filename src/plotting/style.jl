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
