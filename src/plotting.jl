using Plots
using Measures
using LaTeXStrings
using Printf


macro varname(arg)
    return string(arg)
end

function heatmap_myway(x_axis, y_axis, h_grid;
    xlabel=@varname(x_axis), ylabel=@varname(y_axis), hlabel=@varname(h_grid),
    clims=(minimum(h_grid), maximum(h_grid)), size=(600, 500))

    default(c=:viridis,
        left_margin=2mm, right_margin=5mm, top_margin=2mm, bottom_margin=2mm,
        clims=clims,
        size=size)

    return heatmap(x_axis, y_axis, h_grid,
        xlabel=xlabel, ylabel=ylabel,
        title=(hlabel * " against " * xlabel * " and " * ylabel))
end

function plot_myway(x_axis, y_grids;
    xlabel=@varname(x_axis), ylabel="Y_thing", size=(600, 500))

    default(left_margin=2mm, right_margin=5mm, top_margin=2mm, bottom_margin=2mm,
        size=size)

    return plot(x_axis, y_grids,
        xlabel=xlabel, ylabel=ylabel,
        title=(ylabel * " against " * xlabel))
end