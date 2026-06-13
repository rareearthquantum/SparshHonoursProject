using Plots
using Measures
using LaTeXStrings
using Printf


function plot_a_intensity_vs_t_z(t_grid, z_grid, a_grid; clims=(0.0, 1.0), size=(600, 400))
    heatmap(t_grid, z_grid, abs2.(a_grid[:, :, end÷2]),
        c=:viridis,
        xlabel="t", ylabel="z", title="\nElectric field against t and z",
        left_margin=5mm, right_margin=10mm, top_margin=5mm, bottom_margin=5mm,
        clims=clims,
        size=size)
end

function plot_a_intensity_vs_y_z(y_grid, z_grid, a_grid; clims=(0.0, 1.0), size=(600, 400))
    heatmap(y_grid, z_grid, abs2.(a_grid[:, end÷2, :]),
        c=:viridis,
        xlabel="y", ylabel="z", title="\nElectric field intensity against y and z",
        left_margin=5mm, right_margin=10mm, top_margin=5mm, bottom_margin=5mm,
        clims=clims,
        size=size)
end

function plotting_a_intensity(t_grid, y_grid, z_grid, a_grid; clims=(0.0, 1.0), size=(600, 600))
    default(c=:viridis,
        left_margin=2mm, right_margin=5mm, top_margin=2mm, bottom_margin=2mm,
        clims=clims,
        size=size)
    ht = heatmap(t_grid, z_grid, abs2.(a_grid[:, :, end÷2]), xlabel="t", ylabel="z", 
    title=@sprintf("Envelope intensity against t and z at y=%.2g",y_grid[end÷2]))
    hy = heatmap(y_grid, z_grid, abs2.(a_grid[:, end÷2, :]), xlabel="y", ylabel="z", 
    title=@sprintf("Envelope intensity against y and z at t=%.2g",t_grid[end÷2]))
    heatmap(ht, hy, layout=(2, 1))
end