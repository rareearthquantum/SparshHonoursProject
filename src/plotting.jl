function plot_a_intensity_vs_t_z(y_grid, z_grid, a_grid; clims=(0.0, 1.0), size=:default)
    heatmap(t_grid, z_grid, abs2.(a_grid[:, :, end÷2]),
        c=:viridis,
        xlabel="t", ylabel="z", title="\nElectric field against t and z",
        left_margin=5mm, right_margin=10mm, top_margin=5mm, bottom_margin=5mm, clims=clims, size=size)
end

function plot_a_intensity_vs_y_z(y_grid, z_grid, a_grid; clims=(0.0, 1.0), size=:default)
    heatmap(y_grid, z_grid, abs2.(a_grid[:, end÷2, :]), c=:viridis,
        xlabel="y", ylabel="z", title="\nElectric field intensity against y and z",
        left_margin=5mm, right_margin=10mm, top_margin=5mm, bottom_margin=5mm,
        clims=clims, size=size)
end
