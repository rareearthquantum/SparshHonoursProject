"""
    _plot_echo_efficiency(result, field_values)

Return echo energy relative to input-face pulse energy versus z, using
the time-by-z `field_values` matrix. The first two temporal pulses define
the input and retrieval centres; the echo is expected at twice their
separation after the input. Requires two pulses, nonzero input energy,
and a time grid covering the integration windows.
"""
function _plot_echo_efficiency(result, field_values)
    z_vec = result.z_vec
    t_vec = result.t_vec

    input = result.cfg.pulses[1][1]
    retrieval = result.cfg.pulses[2][1]
    tau = retrieval.center - input.center
    echo_time = input.center + 2tau

    input_peaks = zeros(length(z_vec))
    echo_peaks = zeros(length(z_vec))

    input_index = searchsortedfirst(t_vec, input.center)
    echo_index = searchsortedfirst(t_vec, echo_time)
    width_index = searchsortedfirst(t_vec, input.width)

    inputrange = clamp.((input_index-3width_index):(input_index+3width_index), 1, length(t_vec))

    for i in 1:length(z_vec)
        input_peaks[i] = sum(field_values[inputrange, i])
        echo_peaks[i] = sum(field_values[(echo_index-2width_index):end, i])
    end

    efficiency = echo_peaks ./ input_peaks[1]

    return plot(z_vec, efficiency; title="Echo efficiency", xlabel="\$z\$",
        ylabel="\$E_{echo}/E_{input,max}\$", size=PLOT_PANEL_SIZE)
end

"""
    plot_2sd_echo_efficiency(result)

Return echo efficiency versus z from intensity summed over transverse
positions. Requires the two-pulse sequence used by
[`_plot_echo_efficiency`](@ref).
"""
function plot_2sd_echo_efficiency(result)
    field_values = _field_2d_summed_over_y(result, abs2)
    return _plot_echo_efficiency(result, field_values)
end

"""
    plot_1sd_echo_efficiency(result)

Return echo efficiency versus z for a field with a singleton transverse
axis, using the two-pulse windows in [`_plot_echo_efficiency`](@ref).
"""
function plot_1sd_echo_efficiency(result)
    field_values = _field_1d(result, abs2)
    return _plot_echo_efficiency(result, field_values)
end
