using Plots

function heatmap_matrix(A)
    return transpose(A)
end

function plot_propagation(result; plot_size=(1500, 1000))
    time_vec = result.time_vec
    z_vec = result.z_vec
    Omega = result.Omega
    P = result.P

    p_limit = maximum(abs, vcat(vec(real.(P)), vec(imag.(P))))
    omega_limit = maximum(abs, vcat(vec(real.(Omega)), vec(imag.(Omega))))
    p_clims = (-p_limit, p_limit)
    omega_clims = (-omega_limit, omega_limit)

    omega_abs2_max = maximum(abs2, Omega)
    omega_abs2_ylim = (0, omega_abs2_max == 0 ? 1 : 1.05omega_abs2_max)
    z_indices = (1, cld(length(z_vec), 2), length(z_vec))

    heatmap_style = (
        framestyle=:box,
        legend=false,
        colorbar=false,
        grid=false,
        tickfontsize=9,
        guidefontsize=10,
        titlefontsize=11,
        titlefontcolor=RGB(0.15, 0.18, 0.20),
        foreground_color_axis=RGB(0.35, 0.38, 0.40),
        background_color=:white,
        margin=2Plots.mm,
    )
    line_style = (
        framestyle=:box,
        legend=false,
        grid=:y,
        gridalpha=0.16,
        minorgrid=false,
        tickfontsize=9,
        guidefontsize=10,
        titlefontsize=11,
        titlefontcolor=RGB(0.15, 0.18, 0.20),
        foreground_color_axis=RGB(0.35, 0.38, 0.40),
        background_color=:white,
        margin=2Plots.mm,
    )

    heatmap_P_abs2 = heatmap(
        time_vec, z_vec, heatmap_matrix(abs2.(P));
        title="Polarisation intensity  |P|²", xlabel="", ylabel="Position  z",
        c=:viridis, heatmap_style...
    )
    heatmap_P_real = heatmap(
        time_vec, z_vec, heatmap_matrix(real.(P));
        title="Polarisation  Re(P)", xlabel="", ylabel="Position  z",
        c=:viridis, clims=p_clims, heatmap_style...
    )
    heatmap_P_imag = heatmap(
        time_vec, z_vec, heatmap_matrix(imag.(P));
        title="Polarisation  Im(P)", xlabel="Time  t", ylabel="Position  z",
        c=:viridis, clims=p_clims, heatmap_style...
    )

    heatmap_Omega_abs2 = heatmap(
        time_vec, z_vec, heatmap_matrix(abs2.(Omega));
        title="Field intensity  |Ω|²", xlabel="", ylabel="",
        c=:viridis, heatmap_style...
    )
    heatmap_Omega_real = heatmap(
        time_vec, z_vec, heatmap_matrix(real.(Omega));
        title="Field  Re(Ω)", xlabel="", ylabel="",
        c=:viridis, clims=omega_clims, heatmap_style...
    )
    heatmap_Omega_imag = heatmap(
        time_vec, z_vec, heatmap_matrix(imag.(Omega));
        title="Field  Im(Ω)", xlabel="Time  t", ylabel="",
        c=:viridis, clims=omega_clims, heatmap_style...
    )

    omega_input = plot(
        time_vec, abs2.(Omega[:, z_indices[1]]);
        xlabel="", ylabel="", title="Input face",
        ylim=omega_abs2_ylim, linewidth=2.2, color=RGB(0.08, 0.36, 0.48), line_style...
    )
    omega_center = plot(
        time_vec, abs2.(Omega[:, z_indices[2]]);
        xlabel="", ylabel="Field intensity  |Ω|²", title="Centre of medium",
        ylim=omega_abs2_ylim, linewidth=2.2, color=RGB(0.08, 0.36, 0.48), line_style...
    )
    omega_output = plot(
        time_vec, abs2.(Omega[:, z_indices[3]]);
        xlabel="Time  t", ylabel="", title="Output face",
        ylim=omega_abs2_ylim, linewidth=2.2, color=RGB(0.08, 0.36, 0.48), line_style...
    )

    return plot(
        heatmap_P_abs2, heatmap_Omega_abs2,
        omega_input,
        heatmap_P_real, heatmap_Omega_real,
        omega_center,
        heatmap_P_imag, heatmap_Omega_imag,
        omega_output;
        size=plot_size,
        layout=grid(3, 3; widths=[0.36, 0.36, 0.28]),
        link=:x,
        left_margin=1Plots.mm,
        right_margin=1Plots.mm,
        top_margin=5Plots.mm,
        bottom_margin=5Plots.mm,
        subplot_spacing=1Plots.mm,
        background_color=:white,
    )
end

function plot_propagation_compact(result; plot_size=(1540, 900))
    time_vec = result.time_vec
    z_vec = result.z_vec
    Omega = result.Omega
    P = result.P

    p_limit = maximum(abs, vcat(vec(real.(P)), vec(imag.(P))))
    omega_limit = maximum(abs, vcat(vec(real.(Omega)), vec(imag.(Omega))))
    p_clims = (-p_limit, p_limit)
    omega_clims = (-omega_limit, omega_limit)

    omega_abs2_ylim = extrema(abs2.(Omega))
    p_abs2_clims = extrema(abs2.(P))
    omega_abs2_clims = omega_abs2_ylim
    z_indices = (1, cld(length(z_vec), 2), length(z_vec))
    viridis_line_color = palette(:viridis, 3)[2]

    p_abs2_ticks = range(p_abs2_clims...; length=4)
    p_component_ticks = range(p_clims...; length=5)
    omega_abs2_ticks = range(omega_abs2_clims...; length=4)
    omega_component_ticks = range(omega_clims...; length=5)

    heatmap_P_abs2 = heatmap(
        time_vec, z_vec, heatmap_matrix(abs2.(P));
        title="Polarisation intensity",
        xlabel="",
        ylabel="z",
        c=:viridis,
        clims=p_abs2_clims,
        colorbar_ticks=p_abs2_ticks,
        colorbar_tickfontsize=7,
    )
    heatmap_P_imag = heatmap(
        time_vec, z_vec, heatmap_matrix(imag.(P));
        title="Polarisation imaginary part",
        xlabel="t",
        ylabel="z",
        c=:viridis,
        clims=p_clims,
        colorbar_ticks=p_component_ticks,
        colorbar_tickfontsize=7,
    )
    heatmap_Omega_abs2 = heatmap(
        time_vec, z_vec, heatmap_matrix(abs2.(Omega));
        title="Field intensity",
        xlabel="",
        ylabel="",
        c=:viridis,
        clims=omega_abs2_clims,
        colorbar_ticks=omega_abs2_ticks,
        colorbar_tickfontsize=7,
    )
    heatmap_Omega_real = heatmap(
        time_vec, z_vec, heatmap_matrix(real.(Omega));
        title="Field real part",
        xlabel="t",
        ylabel="",
        c=:viridis,
        clims=omega_clims,
        colorbar_ticks=omega_component_ticks,
        colorbar_tickfontsize=7,
    )

    omega_input = plot(
        time_vec, abs2.(Omega[:, z_indices[1]]);
        title="Input face",
        xlabel="",
        ylabel="",
        legend=false,
        color=viridis_line_color,
        ylim=omega_abs2_ylim,
    )
    omega_output = plot(
        time_vec, abs2.(Omega[:, z_indices[3]]);
        title="Output face",
        xlabel="t",
        ylabel="",
        legend=false,
        color=viridis_line_color,
        ylim=omega_abs2_ylim,
    )

    return plot(
        heatmap_P_abs2, heatmap_Omega_abs2, omega_input,
        heatmap_P_imag, heatmap_Omega_real, omega_output;
        size=plot_size,
        layout=(2, 3),
        link=:x,
        left_margin=2Plots.mm,
        right_margin=2Plots.mm,
        top_margin=3Plots.mm,
        bottom_margin=3Plots.mm,
    )
end

function animate_field_and_polarisation(result; filename="field_and_polarisation.gif", fps=30)
    time_vec = result.time_vec
    z_vec = result.z_vec
    Omega = result.Omega
    P = result.P

    omega_ylim = extrema(abs2.(Omega))
    p_ylim = extrema(abs2.(P))

    anim = @animate for i in eachindex(time_vec)
        p1 = plot(
            z_vec,
            abs2.(Omega[i, :]);
            ylabel="abs2(Ω)",
            ylim=omega_ylim,
            legend=false,
            title="t = $(round(time_vec[i], digits=3))",
        )

        p2 = plot(
            z_vec,
            abs2.(P[i, :]);
            xlabel="z",
            ylabel="abs2(P)",
            ylim=p_ylim,
            legend=false,
        )

        plot(p1, p2; layout=(2, 1), size=(700, 600))
    end

    return gif(anim, filename; fps=fps)
end
