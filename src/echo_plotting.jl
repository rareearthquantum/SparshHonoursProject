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

function finite_maximum_abs(values)
    magnitude = 0.0
    for value in values
        isfinite(value) || continue
        magnitude = max(magnitude, abs(value))
    end
    return magnitude
end

function display_scale(values, units; default=(1.0, last(first(units))))
    magnitude = finite_maximum_abs(values)
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

function compact_tick(value)
    isfinite(value) || return string(value)
    iszero(value) && return "0"
    magnitude = abs(value)
    if 1e-2 <= magnitude < 1e4
        return string(round(value; sigdigits=4))
    end

    exponent = floor(Int, log10(magnitude))
    mantissa = round(value / 10.0^exponent; sigdigits=3)
    return "$(mantissa)e$(exponent)"
end

function engineering_scale(values)
    magnitude = finite_maximum_abs(values)
    (iszero(magnitude) || 1e-2 <= magnitude < 1e4) &&
        return (factor=1.0, exponent=0)

    exponent = 3floor(Int, log10(magnitude) / 3)
    return (factor=10.0^exponent, exponent)
end

const SUPERSCRIPT_DIGITS = Dict(
    '-' => '⁻', '0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴',
    '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹',
)
superscript(number::Integer) = String([SUPERSCRIPT_DIGITS[character] for character in string(number)])

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

"""
    plot_2d_echo_snapshots(result; snapshot_times=nothing, plot_size=nothing)

Plot transverse field-intensity snapshots from a 2D propagation result. Each row
contains the field throughout the medium as a `y`-versus-`z` heatmap, followed
by transverse profiles at the input and output faces.

By default, a one-pulse configuration produces one row. A configuration with
two or more pulses produces rows at the first pulse, last pulse, and expected
echo (`2t_last - t_first`). Pass `snapshot_times` to choose times explicitly.
"""
function plot_2d_echo_snapshots(result; snapshot_times=nothing, plot_size=nothing)
    ndims(result.Omega) == 3 || throw(ArgumentError(
        "plot_2d_echo_snapshots requires Omega indexed as (time, z, y)"))

    y_vec = hasproperty(result, :y_vec) ? result.y_vec : make_y_grid(result.cfg)
    size(result.Omega) == (length(result.time_vec), length(result.z_vec), length(y_vec)) ||
        throw(DimensionMismatch(
            "Omega dimensions must match the time, z, and y coordinate vectors"))

    automatic_times = isnothing(snapshot_times)
    if automatic_times
        isempty(result.cfg.pulses) && throw(ArgumentError(
            "the configuration has no pulses; pass snapshot_times explicitly"))
        time_pulses = first.(result.cfg.pulses)
        first_time = first(time_pulses).center
        snapshot_times = if length(time_pulses) == 1
            [first_time]
        else
            last_time = last(time_pulses).center
            [first_time, last_time, 2last_time - first_time]
        end
    else
        snapshot_times = collect(snapshot_times)
        isempty(snapshot_times) && throw(ArgumentError("snapshot_times cannot be empty"))
    end

    time_min, time_max = extrema(result.time_vec)
    all(t -> time_min <= t <= time_max, snapshot_times) || throw(ArgumentError(
        "all snapshot times must lie in the simulation window [$time_min, $time_max]"))
    time_indices = [argmin(abs.(result.time_vec .- t)) for t in snapshot_times]

    nonfinite_count = count(value -> !isfinite(value), result.Omega)
    if nonfinite_count > 0
        @warn "Omega contains $nonfinite_count NaN or Inf values; these will appear as gaps in the plot. The propagation result is numerically unstable."
    end

    scales = plot_scales(result)
    y_scale = display_scale(y_vec, LENGTH_DISPLAY_UNITS; default=(1.0, "m"))
    y = y_vec ./ y_scale.factor
    z = result.z_vec ./ scales.position.factor
    intensity = abs2.(result.Omega) ./ scales.field.factor^2
    intensity[.!isfinite.(intensity)] .= NaN
    intensity_scale = engineering_scale(@view intensity[time_indices, :, :])
    intensity ./= intensity_scale.factor
    intensity_units = intensity_scale.exponent == 0 ?
        "($(scales.field.unit))²" :
        "10$(superscript(intensity_scale.exponent)) ($(scales.field.unit))²"
    intensity_label = "|Ω|² ($intensity_units)"
    y_label = axis_label("Transverse position, y", y_scale)
    z_label = axis_label("Propagation distance, z", scales.position)

    row_names = if !automatic_times
        ["Snapshot $i" for i in eachindex(snapshot_times)]
    elseif length(snapshot_times) == 1
        ["Pulse"]
    else
        ["First pulse", "Last pulse", "Expected echo"]
    end

    panel_style = (
        framestyle=:box,
        tickfontsize=9,
        guidefontsize=10,
        titlefontsize=11,
        margin=2Plots.mm,
    )
    line_style = (
        panel_style...,
        linewidth=2.2,
        legend=false,
        grid=:y,
        gridalpha=0.16,
        xformatter=compact_tick,
        yformatter=compact_tick,
    )

    panels = Plots.Plot[]
    for (row, time_index) in enumerate(time_indices)
        shown_time = result.time_vec[time_index] / scales.time.factor
        time_text = "$(row_names[row]) — t = $(compact_tick(shown_time)) $(scales.time.unit)"
        bottom_row = row == length(time_indices)

        push!(panels, heatmap(y, z, @view(intensity[time_index, :, :]);
            title="$time_text\nField intensity in medium",
            xlabel=bottom_row ? y_label : "",
            ylabel=z_label,
            xformatter=compact_tick,
            yformatter=compact_tick,
            c=:viridis,
            clims=extrema(intensity),
            left_margin=6Plots.mm,
            right_margin=5Plots.mm,
            panel_style...,
        ))
        push!(panels, plot(y, @view(intensity[time_index, firstindex(z), :]);
            title=row == 1 ? "Input face" : "",
            xlabel=bottom_row ? y_label : "",
            ylabel=intensity_label,
            color=palette(:viridis, 3)[2],
            clims=extrema(intensity),
            line_style...,
        ))
        push!(panels, plot(y, @view(intensity[time_index, lastindex(z), :]);
            title=row == 1 ? "Output face" : "",
            xlabel=bottom_row ? y_label : "",
            ylabel=intensity_label,
            color=palette(:viridis, 3)[3],
            clims=clims=extrema(intensity),
            line_style...,
        ))
    end

    nrows = length(time_indices)
    isnothing(plot_size) && (plot_size = (1500, max(420, 340nrows)))
    return plot(panels...;
        layout=(nrows, 3),
        size=plot_size,
        left_margin=6Plots.mm,
        right_margin=5Plots.mm,
        top_margin=6Plots.mm,
        bottom_margin=5Plots.mm,
        background_color=:white,
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

function plot_echo_transmission_and_efficiency_vs_optical_depth(result)
    echo_time = 2*result.cfg.pulses[2].center - result.cfg.pulses[1].center
    echo_time_index = time_index_grabber(echo_time, result.time_vec)

    input_pulse_time = result.cfg.pulses[1].center
    retrieval_pulse_time = result.cfg.pulses[2].center
    echo_time = 2*retrieval_pulse_time - input_pulse_time
    halfway_retrieval_echo_time = (echo_time + retrieval_pulse_time) / 2
    tieme = time_index_grabber(halfway_retrieval_echo_time, result.time_vec)

    echo_intensity_maxes = [maximum(abs2.(result.Omega[tieme:end,i])) for i in eachindex(result.z_vec)]

    input_pulse_max = result.cfg.pulses[1].area
    efficiencies = echo_intensity_maxes./input_pulse_max

    fig = plot(result.z_vec, echo_intensity_maxes, label="Transmission")
    plot!(result.z_vec, efficiencies, label="Efficiency")

    return fig
end


function plot_propagation_super_compact(result; plot_size=(1540, 900))
    scales = plot_scales(result)
    time = result.time_vec ./ scales.time.factor
    position = result.z_vec ./ scales.position.factor
    time_label = axis_label("Time, t", scales.time)
    position_label = axis_label("Position, z", scales.position)
    field_unit = scales.field.unit

    P_abs2 = abs2.(result.P)
    Omega_abs2 = abs2.(result.Omega) ./ scales.field.factor^2

    p_abs2_clims = extrema(P_abs2)
    omega_abs2_clims = extrema(Omega_abs2)
    z_indices = (1, length(position))
    heatmap_style = (colorbar_tickfontsize=7,)

    panels = (
        heatmap_panel(time, position, Omega_abs2;
            title="Field intensity |Ω|² ($field_unit)²",
            clims=omega_abs2_clims,
            colorbar_ticks=range(omega_abs2_clims...; length=4),
            heatmap_style...),
        heatmap_panel(time, position, P_abs2;
            title="Polarisation intensity |P|²", ylabel=position_label,
            clims=p_abs2_clims, colorbar_ticks=range(p_abs2_clims...; length=4),
            heatmap_style...),
        field_panel(time, result.Omega, z_indices[1], scales.field.factor;
            title="Input face", ylim=omega_abs2_clims),
        field_panel(time, result.Omega, z_indices[2], scales.field.factor;
            title="Output face", xlabel=time_label, ylim=omega_abs2_clims),
    )

    line_plots = plot(panels[end-1], panels[end],layout=(2,1))
    heatmaps = plot(panels[begin], panels[begin+1],layout=(2,1))

    return plot(heatmaps,line_plots;
        size=plot_size,
        layout=(1, 2),
        link=:x,
        left_margin=2Plots.mm,
        right_margin=2Plots.mm,
        top_margin=3Plots.mm,
        bottom_margin=3Plots.mm,
    )
end