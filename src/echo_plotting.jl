using Plots

const TIME_DISPLAY_UNITS = (
    (1.0, "s"), (1e-3, "ms"), (1e-6, "μs"), (1e-9, "ns"), (1e-12, "ps"),
)
const LENGTH_DISPLAY_UNITS = (
    (1e3, "km"), (1.0, "m"), (1e-2, "cm"),
    (1e-3, "mm"), (1e-6, "μm"), (1e-9, "nm"),
)
const RATE_DISPLAY_UNITS = (
    (1e12, "ps⁻¹"), (1e9, "ns⁻¹"), (1e6, "μs⁻¹"),
    (1e3, "ms⁻¹"), (1.0, "s⁻¹"),
)

function display_scale(values, units; default=(1.0, last(first(units))))
    magnitude = maximum(abs, values)
    iszero(magnitude) && return (factor=default[1], unit=default[2])

    for (factor, unit) in units
        magnitude >= factor && return (; factor, unit)
    end

    factor, unit = last(units)
    return (; factor, unit)
end

function plot_scales(result)
    time = display_scale(result.time_vec, TIME_DISPLAY_UNITS; default=(1.0, "s"))
    position = display_scale(result.z_vec, LENGTH_DISPLAY_UNITS; default=(1.0, "m"))
    field = display_scale(result.Omega, RATE_DISPLAY_UNITS; default=(1.0, "s⁻¹"))
    return (; time, position, field)
end

axis_label(name, scale) = "$name ($(scale.unit))"

function symmetric_limits(arrays...)
    limit = maximum(maximum(abs, array) for array in arrays)
    iszero(limit) && (limit = 1)
    return (-limit, limit)
end

function intensity_ylim(values)
    upper = maximum(values)
    return (0, iszero(upper) ? 1 : 1.05upper)
end

function heatmap_panel(x, y, values; kwargs...)
    return heatmap(x, y, transpose(values); c=:viridis, kwargs...)
end

function field_panel(time, Omega, z_index, field_factor; kwargs...)
    intensity = @views abs2.(Omega[:, z_index]) ./ field_factor^2
    return plot(time, intensity; color=palette(:viridis, 3)[2], legend=false, kwargs...)
end

function plot_propagation(result; plot_size=(1500, 1000))
    scales = plot_scales(result)
    time = result.time_vec ./ scales.time.factor
    position = result.z_vec ./ scales.position.factor
    time_label = axis_label("Time, t", scales.time)
    position_label = axis_label("Position, z", scales.position)
    field_unit = scales.field.unit

    P_abs2 = abs2.(result.P)
    P_real = real.(result.P)
    P_imag = imag.(result.P)
    Omega_abs2 = abs2.(result.Omega) ./ scales.field.factor^2
    Omega_real = real.(result.Omega) ./ scales.field.factor
    Omega_imag = imag.(result.Omega) ./ scales.field.factor

    p_clims = symmetric_limits(P_real, P_imag)
    omega_clims = symmetric_limits(Omega_real, Omega_imag)
    omega_ylim = intensity_ylim(Omega_abs2)
    z_indices = (1, cld(length(position), 2), length(position))

    panel_style = (
        framestyle=:box,
        tickfontsize=9,
        guidefontsize=10,
        titlefontsize=11,
        titlefontcolor=RGB(0.15, 0.18, 0.20),
        foreground_color_axis=RGB(0.35, 0.38, 0.40),
        background_color=:white,
        margin=2Plots.mm,
    )
    heatmap_style = (; panel_style..., colorbar=false, grid=false)
    line_style = (; panel_style..., grid=:y, gridalpha=0.16, minorgrid=false,
        linewidth=2.2, ylim=omega_ylim)

    panels = (
        heatmap_panel(time, position, P_abs2;
            title="Polarisation intensity |P|²", ylabel=position_label,
            heatmap_style...),
        heatmap_panel(time, position, Omega_abs2;
            title="Field intensity |Ω|² ($field_unit)²", heatmap_style...),
        field_panel(time, result.Omega, z_indices[1], scales.field.factor;
            title="Input face", line_style...),
        heatmap_panel(time, position, P_real;
            title="Polarisation Re(P)", ylabel=position_label, clims=p_clims,
            heatmap_style...),
        heatmap_panel(time, position, Omega_real;
            title="Field Re(Ω) ($field_unit)", clims=omega_clims,
            heatmap_style...),
        field_panel(time, result.Omega, z_indices[2], scales.field.factor;
            title="Centre of medium", ylabel="Field intensity ($field_unit)²",
            line_style...),
        heatmap_panel(time, position, P_imag;
            title="Polarisation Im(P)", xlabel=time_label, ylabel=position_label,
            clims=p_clims, heatmap_style...),
        heatmap_panel(time, position, Omega_imag;
            title="Field Im(Ω) ($field_unit)", xlabel=time_label,
            clims=omega_clims, heatmap_style...),
        field_panel(time, result.Omega, z_indices[3], scales.field.factor;
            title="Output face", xlabel=time_label, line_style...),
    )

    return plot(panels...;
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
    scales = plot_scales(result)
    time = result.time_vec ./ scales.time.factor
    position = result.z_vec ./ scales.position.factor
    time_label = axis_label("Time, t", scales.time)
    position_label = axis_label("Position, z", scales.position)
    field_unit = scales.field.unit

    P_abs2 = abs2.(result.P)
    P_imag = imag.(result.P)
    Omega_abs2 = abs2.(result.Omega) ./ scales.field.factor^2
    Omega_real = real.(result.Omega) ./ scales.field.factor

    p_clims = symmetric_limits(P_imag)
    omega_clims = symmetric_limits(Omega_real)
    p_abs2_clims = extrema(P_abs2)
    omega_abs2_clims = extrema(Omega_abs2)
    z_indices = (1, length(position))
    heatmap_style = (colorbar_tickfontsize=7,)

    panels = (
        heatmap_panel(time, position, P_abs2;
            title="Polarisation intensity |P|²", ylabel=position_label,
            clims=p_abs2_clims, colorbar_ticks=range(p_abs2_clims...; length=4),
            heatmap_style...),
        heatmap_panel(time, position, Omega_abs2;
            title="Field intensity |Ω|² ($field_unit)²",
            clims=omega_abs2_clims,
            colorbar_ticks=range(omega_abs2_clims...; length=4),
            heatmap_style...),
        field_panel(time, result.Omega, z_indices[1], scales.field.factor;
            title="Input face", ylim=omega_abs2_clims),
        heatmap_panel(time, position, P_imag;
            title="Polarisation Im(P)", xlabel=time_label, ylabel=position_label,
            clims=p_clims, colorbar_ticks=range(p_clims...; length=5),
            heatmap_style...),
        heatmap_panel(time, position, Omega_real;
            title="Field Re(Ω) ($field_unit)", xlabel=time_label,
            clims=omega_clims, colorbar_ticks=range(omega_clims...; length=5),
            heatmap_style...),
        field_panel(time, result.Omega, z_indices[2], scales.field.factor;
            title="Output face", xlabel=time_label, ylim=omega_abs2_clims),
    )

    return plot(panels...;
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
    scales = plot_scales(result)
    time = result.time_vec ./ scales.time.factor
    position = result.z_vec ./ scales.position.factor
    position_label = axis_label("Position, z", scales.position)
    field_intensity = abs2.(result.Omega) ./ scales.field.factor^2
    polarisation_intensity = abs2.(result.P)

    anim = @animate for i in eachindex(time)
        field_plot = plot(position, field_intensity[i, :];
            ylabel="|Ω|² ($(scales.field.unit))²",
            ylim=extrema(field_intensity),
            legend=false,
            title="t = $(round(time[i]; sigdigits=4)) $(scales.time.unit)",
        )
        polarisation_plot = plot(position, polarisation_intensity[i, :];
            xlabel=position_label,
            ylabel="|P|²",
            ylim=extrema(polarisation_intensity),
            legend=false,
        )
        plot(field_plot, polarisation_plot; layout=(2, 1), size=(700, 600))
    end

    return gif(anim, filename; fps)
end
