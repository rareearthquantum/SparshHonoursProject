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
