using Plots, Measures
using LsqFit

const PLOT_PANEL_SIZE = (800, 600)
const PLOT_LEFT_MARGIN = 2mm
const PLOT_BOTTOM_MARGIN = 2mm

function _plot_margins(layout)
    layout == (1, 2) && return (5mm, 3mm)
    layout == (1, 3) && return (3mm, 4mm)
    layout == (2, 2) && return (2mm, 1mm)
    return (PLOT_LEFT_MARGIN, PLOT_BOTTOM_MARGIN)
end

function _normalise_by_max!(values)
    values ./= maximum(values)
    return values
end

function _unit_scale_ylims(values; padding_fraction=0.05)
    lower = min(0.0, minimum(values))
    upper = max(1.0, maximum(values))
    span = upper - lower
    padded_lower = lower < 0 ? lower - padding_fraction * span : 0.0
    padded_upper = upper + padding_fraction * span
    return (padded_lower, padded_upper)
end

_field_2d_summed_over_y(result, operation) =
    dropdims(sum(operation.(result.Omega); dims=3); dims=3)

_field_1d(result, operation) = operation.(dropdims(result.Omega; dims=3))

_slice_colors(count) = cgrad([:blue, :yellow, :red])[range(0, 1, length=count)]

function _slice_label(index, count)
    index == 1 && return "first"
    index == count && return "last"
    index == count ÷ 2 + 1 && return "middle"
    return false
end

function _add_panel_titles!(plots)
    for (index, panel) in enumerate(plots)
        panel_label = "($(Char(Int('a') + index - 1)))"
        title!(panel, panel_label; titlelocation=:left)
    end
    return plots
end

function _compose_plots(plots...; layout)
    _add_panel_titles!(plots)
    rows, columns = layout
    height = round(Int, PLOT_PANEL_SIZE[2] * rows / columns)
    size = (PLOT_PANEL_SIZE[1], height)
    left_margin, bottom_margin = _plot_margins(layout)
    return plot(plots...; layout=layout, size=size,
        left_margin=left_margin, bottom_margin=bottom_margin)
end

function _plot_tz_heatmap(time, z, field; title)
    return heatmap(time, z, transpose(field);
        c=:viridis, xlabel="\$t\$", ylabel="\$z\$", title=title,
        cbar=false, size=PLOT_PANEL_SIZE)
end

function _plot_xlims(x, field, zoom)
    zoom || return extrema(x)
    startindex, endindex = zoomingin(field)
    return (x[startindex], x[endindex])
end

function _plot_temporal_slices(field, time, indices; title, zoom, ylabel, size=PLOT_PANEL_SIZE)
    xlims = _plot_xlims(time, field, zoom)
    colors = _slice_colors(length(indices))
    fig = plot(; xlabel="\$t\$", ylabel=ylabel, title=title, xlims=xlims, size=size)

    for (i, field_index) in enumerate(indices)
        plot!(time, field[:, field_index]; color=colors[i], label=_slice_label(i, length(indices)))
    end

    return fig
end

_gaussian_fit(x, p) = @. p[1] * exp(-0.5 * ((x-p[2])/p[3])^2) + p[4]
_sech_squared_fit(x, p) = @. p[1] / cosh((x-p[2])/p[3])^2 + p[4]

function _fit_initial_parameters(x, y)
    return [maximum(y), x[argmax(y)], (maximum(x)-minimum(x))/10, 0.0]
end

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

#--------------------------

function plot_2d_tz_ysum_heatmap(result; operation=abs2, title=:default)
    Omega = _field_2d_summed_over_y(result, operation)
    _normalise_by_max!(Omega)
    title = (title==:default) ? "sum $(nameof(operation)) Omega over y" : title
    return _plot_tz_heatmap(result.time_vec, result.z_vec, Omega; title=title)
end

function plot_2d_z_ysum_line(result; operation=abs2, title=:default,
        size=PLOT_PANEL_SIZE, ylabel="", label=false, unit_scale=false)
    Omega = _field_2d_summed_over_y(result, operation)
    z_vec = result.z_vec

    Omega_z = vec(sum(Omega; dims=1))
    Omega_z ./= Omega_z[begin]

    x_axis = z_vec
    xlims = extrema(x_axis)

    title = (title==:default) ? "sum $(nameof(operation)) Omega over y and t" : title

    fig = plot(x_axis, Omega_z; xlabel="\$z\$", ylabel=ylabel,
        title=title, label=label, xlims=xlims, size=size)
    unit_scale && ylims!(fig, _unit_scale_ylims(Omega_z))
    return fig
end

function plot_2d_superimposed_z_ysum_line(result; operation=abs2, nslices=10, title=:default, zoom=true, ylabel="")
    Omega = _field_2d_summed_over_y(result, operation)
    _normalise_by_max!(Omega)

    title = (title==:default) ? "sum $(nameof(operation)) Omega over y" : title
    delta = floor(Int, length(result.z_vec) / nslices)
    indices = 1 .+ (0:nslices-1) .* delta
    return _plot_temporal_slices(Omega, result.time_vec, indices;
        title=title, zoom=zoom, ylabel=ylabel)
end


function plot_2d_z_yslices_line(result; operation=abs2, nslices=10, title=:default, ylabel="")
    Omega = operation.(result.Omega[:, :, (end÷2+1):end])
    z_vec = result.z_vec
    title = (title==:default) ? "sum $(nameof(operation)) Omega and t" : title
    Omega_z = dropdims(sum(Omega; dims=1); dims=1)
    _normalise_by_max!(Omega_z)
    x_axis = z_vec
    xlims = extrema(x_axis)

    colors = _slice_colors(nslices)
    fig = plot(; xlabel="\$z\$", ylabel=ylabel, title=title, xlims=xlims,
        legend=false, size=PLOT_PANEL_SIZE)

    delta = floor(Int, length(result.y_vec) / 2nslices)


    for i in 1:nslices
        plot!(x_axis, Omega_z[:, i*delta-(delta-1)], c=colors[i])
    end

    return fig
end

function plot_2d_zy_pulseprofile_heatmap(result; operation=abs2, title=:default)
    Omega = operation.(result.Omega)
    _normalise_by_max!(Omega)
    y_vec, z_vec = result.y_vec, result.z_vec

    tindex = argmax(Omega[:, end÷2+1, end÷2+1])
    Omega_tslice = Omega[tindex, :, :]
    _normalise_by_max!(Omega_tslice)

    title = (title==:default) ? "$(nameof(operation)) Omega time slice" : title

    return heatmap(y_vec, z_vec, Omega_tslice; c=:viridis,
        xlabel="\$y\$", ylabel="\$z\$", title=title,
        cbar=false, size=PLOT_PANEL_SIZE)
end


function plot_2d_superimposed_lines_against_t(result; operation=abs2, nslices=10, title="", zoom=false, ylabel="")
    Omega = _field_2d_summed_over_y(result, operation)
    _normalise_by_max!(Omega)

    title = (title==:default) ? "$(nameof(operation)) Omega" : title
    delta = floor(Int, length(result.z_vec) / nslices)
    indices = 1 .+ (0:nslices-1) .* delta
    return _plot_temporal_slices(Omega, result.time_vec, indices;
        title=title, zoom=zoom, ylabel=ylabel)
end

function _plot_echo_efficiency(result, Omega)
    z_vec = result.z_vec
    t_vec = result.time_vec

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
        input_peaks[i] = sum(Omega[inputrange, i])
        echo_peaks[i] = sum(Omega[(echo_index-2width_index):end, i])
    end

    efficiency = echo_peaks ./ input_peaks[1]

    return plot(z_vec, efficiency; title="Echo efficiency", xlabel="\$z\$",
        ylabel="\$E_{echo}/E_{input,max}\$", size=PLOT_PANEL_SIZE)
end

function plot_2sd_echo_efficiency(result)
    Omega = _field_2d_summed_over_y(result, abs2)
    return _plot_echo_efficiency(result, Omega)
end

function plot_2d_superimposed_lines_against_y(result; operation=abs2, nslices=10, title=:default, size=PLOT_PANEL_SIZE, ylabel="")
    Omega = operation.(result.Omega[end÷2+1, :, :])
    _normalise_by_max!(Omega)

    y_vec = result.y_vec
    z_vec = result.z_vec
    title = (title==:default) ? "$(nameof(operation)) Omega" : title
    x_axis = y_vec
    xlims = extrema(x_axis) #watch it

    colors = _slice_colors(nslices)
    fig = plot(; xlabel="\$y\$", ylabel=ylabel, title=title, xlims=xlims, size=size)

    indices = round.(Int, range(1, length(z_vec), length=nslices))

    for (i, j) in enumerate(indices)
        plot!(x_axis, Omega[j, :], c=colors[i], label=_slice_label(i, nslices))
    end

    return fig
end

function plots_for_jevon_2sd_echo(result)
    heatmap = plot_2d_tz_ysum_heatmap(result; operation=abs2, title="Pulse intensity")
    superimposed_zlines_against_t = plot_2d_superimposed_lines_against_t(result; title="Pulse intensity slices in \$z\$", zoom=true, nslices=10, ylabel="\$I/I_{max}\$")
    echo_efficiency = plot_2sd_echo_efficiency(result)

    return _compose_plots(superimposed_zlines_against_t, echo_efficiency;
        layout=(1, 2))
end

function plots_for_jevon_2sd_SIT(result)
    heatmap = plot_2d_tz_ysum_heatmap(result; operation=abs2, title="Pulse intensity")
    #area = plot_2d_z_ysum_line(result; operation=real, title="Total pulse area")
    pulse_profile_heatmap = plot_2d_zy_pulseprofile_heatmap(result; operation=abs2, title="Pulse intensity profile time slice")
    energy = plot_2d_z_ysum_line(result; operation=abs2,
        title="Total pulse energy", ylabel="\$E/E(0)\$", unit_scale=true)
    superimposed = plot_2d_superimposed_z_ysum_line(result; operation=abs2, nslices=20, title="Pulse intensity slices in \$z\$", zoom=true, ylabel="\$I/I_{max}\$")
    superimposed_area_yslices = plot_2d_z_yslices_line(result; operation=real, nslices=30, title="Pulse area slices in \$y\$ against \$z\$", ylabel="\$A/A_{max}\$")
    superimposed_intensity_yslices = plot_2d_z_yslices_line(result; operation=abs2, nslices=30, title="Pulse intensity slices in \$y\$ against \$z\$", ylabel="\$I/I_{max}\$")

    return _compose_plots(
        heatmap, pulse_profile_heatmap, superimposed_area_yslices,
        superimposed, energy, superimposed_intensity_yslices,
        layout=(2, 3))
end


#-------------------------------------------

function plot_1d_tz_heatmap(result; operation=abs2, title=:default)
    Omega = _field_1d(result, operation)
    _normalise_by_max!(Omega)
    title = (title==:default) ? "$(nameof(operation)) Omega" : title
    return _plot_tz_heatmap(result.time_vec, result.z_vec, Omega; title=title)
end

function plot_1d_superimposed_lines_against_t(result; operation=abs2, nslices=10, title=:default, zoom=false, size=PLOT_PANEL_SIZE, ylabel="")
    Omega = _field_1d(result, operation)
    _normalise_by_max!(Omega)

    title = (title==:default) ? "$(nameof(operation)) Omega" : title
    indices = round.(Int, range(1, length(result.z_vec), length=nslices))
    return _plot_temporal_slices(Omega, result.time_vec, indices;
        title=title, zoom=zoom, ylabel=ylabel, size=size)
end

function plot_1d_superimposed_lines_against_z(result; operation=abs2, nslices=10, title=:default)
    Omega = _field_1d(result, operation)
    _normalise_by_max!(Omega)

    t_vec, z_vec = result.time_vec, result.z_vec
    title = (title==:default) ? "$(nameof(operation)) Omega" : title
    x_axis = z_vec
    xlims = extrema(x_axis) #watch it

    colors = _slice_colors(nslices)
    fig = plot(; xlabel="\$z\$", title=title, xlims=xlims, size=PLOT_PANEL_SIZE)

    delta = floor(Int, length(t_vec) / nslices)


    for i in 1:nslices
        plot!(x_axis, Omega[i*delta-(delta-1), :], c=colors[i], label=_slice_label(i, nslices))
    end

    return fig
end

function plot_1sd_echo_efficiency(result)
    Omega = _field_1d(result, abs2)
    return _plot_echo_efficiency(result, Omega)
end

function plots_for_jevon_1sd_echo(result)
    heatmap = plot_1d_tz_heatmap(result; operation=abs2, title="Pulse intensity")
    superimposed_zlines_against_t = plot_1d_superimposed_lines_against_t(result; title="Pulse intensity slices in \$z\$", zoom=true, nslices=10, ylabel="\$I/I_{max}\$")
    echo_efficiency = plot_1sd_echo_efficiency(result)

    return _compose_plots(heatmap, superimposed_zlines_against_t, echo_efficiency;
        layout=(1, 3))
end

function plots_for_jevon_1sd_SIT(result)
    heatmap = plot_1d_tz_heatmap(result; operation=abs2, title="Pulse intensity")
    superimposed_zlines_against_t = plot_1d_superimposed_lines_against_t(result; title="Pulse intensity slices in \$z\$", zoom=true, nslices=10, ylabel="\$I/I_{max}\$")
    energy = plot_2d_z_ysum_line(result; operation=abs2,
        title="Total pulse energy", ylabel="\$E/E(0)\$", unit_scale=true)
    area = plot_2d_z_ysum_line(result; operation=real,
        title="Pulse area", ylabel="\$A/A(0)\$", unit_scale=true)
    #superimposed_tlines_against_z = plot_1d_superimposed_lines_against_z(result; operation=abs2, nslices=10, title="Pulse intensity slices in \$t\$")

    return _compose_plots(
        heatmap, energy,
        superimposed_zlines_against_t, area,
        layout=(2, 2))
end

#------------------------

function plot_0sd_polarisation(result)
    P = abs2.(result.P)

    t_vec = result.time_vec
    #title = (title==:default) ? "" : title
    fig = plot(; xlabel="\$t\$", ylabel="\$|\\mathcal{P}|^2\$",
        title="Squared magnitude of the polarisation density against time",
        size=PLOT_PANEL_SIZE)

    plot!(t_vec, P; label=false)

    return fig
end

function plot_0sd_pol_complex_parts(result)
    Pu = 2*real(result.P)
    Pv = -2*imag(result.P)

    t_vec = result.time_vec

    fig = plot(; xlabel="\$t\$", ylabel="Bloch component",
        title="Ensemble averaged Bloch components against time",
        size=PLOT_PANEL_SIZE)

    plot!(t_vec, Pu; label="averaged u")
    plot!(t_vec, Pv; label="averaged v")

    return fig
end

function plots_for_jevon_0sd(result)
    pol_abs2 = plot_0sd_polarisation(result)
    pol_real_n_imag = plot_0sd_pol_complex_parts(result)


    return _compose_plots(pol_abs2, pol_real_n_imag; layout=(2, 1))
end

#--------------------------------------------

function plot_fitted_SIT(result; zoom=true)
    Omega = _field_2d_summed_over_y(result, abs2)
    _normalise_by_max!(Omega)
    t_vec = result.time_vec

    # Your simulation data
    x = t_vec
    y1 = Omega[:, begin]
    y2 = Omega[:, end]

    # Initial guesses
    p01 = _fit_initial_parameters(x, y1)
    p02 = _fit_initial_parameters(x, y2)

    # Fit both
    fit_g = curve_fit(_gaussian_fit, x, y1, p01)
    fit_s = curve_fit(_sech_squared_fit, x, y2, p02)

    # Best-fitting parameters
    println(fit_g.param)
    println(fit_s.param)

    gaussian_curve = _gaussian_fit(x, fit_g.param)
    sech_curve = _sech_squared_fit(x, fit_s.param)
    xlims = _plot_xlims(t_vec, Omega, zoom)

    fig = plot(x, y1; title="", ylabel="\$I/I_{max}\$", xlabel="\$t\$", label="Input face data", xlims=xlims, legend=:topleft, st=:scatter, shape=:circle, ms=3, mc=:red, msc=:red)
    plot!(x, y2; label="Output face data", st=:scatter, shape=:circle, ms=3, mc=:blue, msc=:blue)
    plot!(x, gaussian_curve, label="Gaussian fit", st=:line, lc=:red, ls=:solid)
    plot!(x, sech_curve, label="Sech^2 fit", st=:line, lc=:blue, ls=:solid)

    difffig = plot(x, y1 - gaussian_curve; label="Input face data-fit difference", xlims=xlims, ylabel="\$I/I_{max}\$", xlabel="\$t\$", legend=:topleft)
    plot!(x, y2 - sech_curve; label="Output face data-fit difference")

    return _compose_plots(fig, difffig; layout=(1, 2))
end

function plot_fitted_attenuation(result)
    Omega = _field_2d_summed_over_y(result, abs2)
    _normalise_by_max!(Omega)
    z_vec = result.z_vec

    # Your simulation data
    x = z_vec
    y = vec(sum(Omega[:, :], dims=1))
    _normalise_by_max!(y)

    z0 = x[begin]
    exponentialfit(x, p) = @. p[1] * exp(-(x-z0)/p[2]) + p[3]

    # Initial guesses
    p0 = [maximum(y), (maximum(x)-minimum(x))/10, 0.0]

    # Fit both
    fit_e = curve_fit(exponentialfit, x, y, p0)
    expfit = exponentialfit(x, fit_e.param)


    fig = plot(x, y; title="", xlabel ="\$z\$", ylabel="\$E/E_{max}\$", label="Data", st=:scatter, shape=:circle, ms=3, mc=:red, msc=:red, lc=:red)
    plot!(x, expfit, label="Exponential fit", st=:line, lc=:red, ls=:solid)

    difffig = plot(x, y - expfit; label="Data-fit difference", xlabel="\$z\$", ylabel="\$E/E_{max}\$")

    return _compose_plots(fig, difffig; layout=(1, 2))
end


function plot_fitted_diffraction(result)
    Omega = abs2.(result.Omega[end÷2+1, :, :]) #halfway in time
    _normalise_by_max!(Omega)
    y_vec = result.y_vec

    # Your simulation data
    x = y_vec
    y1 = Omega[begin, :]
    y2 = Omega[end, :]

    # Initial guesses
    p01 = _fit_initial_parameters(x, y1)
    p02 = _fit_initial_parameters(x, y2)

    # Fit both
    fit_g1 = curve_fit(_gaussian_fit, x, y1, p01)
    fit_g2 = curve_fit(_gaussian_fit, x, y2, p02)

    # Best-fitting parameters
    println(fit_g1.param)
    println(fit_g2.param)

    input_curve = _gaussian_fit(x, fit_g1.param)
    output_curve = _gaussian_fit(x, fit_g2.param)

    fig = plot(x, y1; title="", ylabel="\$I/I_{max}\$", label="Input face data", xlabel="\$y\$", st=:scatter, shape=:circle, ms=3, mc=:red, msc=:red, lc=:red)
    plot!(x, y2; label="Output face data", st=:scatter, shape=:circle, ms=3, mc=:blue, msc=:blue, lc=:blue)
    plot!(x, input_curve, label="Fit for input face", st=:line, lc=:red, ls=:solid)
    plot!(x, output_curve, label="Fit for output face", st=:line, lc=:blue, ls=:solid)

    difffig = plot(x, y1 - input_curve; label="Input face data-fit difference", ylabel="\$I/I_{max}\$", xlabel="\$y\$")
    plot!(x, y2 - output_curve; label="Output face data-fit difference")

    return _compose_plots(fig, difffig; layout=(1, 2))
end
