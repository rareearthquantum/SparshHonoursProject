using Plots, Measures


function animate_field_2d(result; filename="echo_2d.gif", fps=30, operation=abs2)
    field_intensity = operation.(result.Omega)
    z_vec = result.z_vec
    y_vec = result.y_vec
    t_vec = result.time_vec

    anim = @animate for i in eachindex(t_vec)
        heatmap(y_vec, z_vec, field_intensity[i, :, :];
            clabel="|Ω|²",
            ylabel="z",
            xlabel="y",
            clims=extrema(field_intensity),
            legend=false,
            title="Rabi frequency intensity at t = $(round(t_vec[i]; sigdigits=4))",
            c=:viridis
        )
    end

    return anim
end

#---------------------------

function zoomingin(vector;threshold=0.001)

    first_above = lastindex(vector[:,1])
    last_above = firstindex(vector[:,1])

    for i in 1:length(vector[1,:])
        viewvec = @view vector[:,i]
        temp = findfirst(>=(threshold), viewvec)
        if temp !== nothing && temp < first_above
            first_above = temp
        end

        temp = findlast(>=(threshold), viewvec)
        if temp !== nothing && temp > last_above
            last_above = temp
        end
    end

    startindex = max(first_above - 1, firstindex(vector[:,1]))
    endindex = min(last_above + 1, lastindex(vector[:,1]))
    
    return (startindex,endindex)
end

#--------------------------

function plot_2d_tz_ysum_heatmap(result; operation=abs2, title=:default)
    Omega = operation.(dropdims(sum(result.Omega; dims=3); dims=3))
    Omega ./= maximum(Omega)
    t_vec, z_vec = result.time_vec, result.z_vec
    title = (title==:default) ? "sum $(nameof(operation)) Omega over y" : title
    return heatmap(t_vec, z_vec, transpose(Omega[:, :]);
        c=:viridis, xlabel="\$t\$", ylabel="\$z\$", title=title, cbar=false)
end

function plot_2d_z_ysum_line(result; operation=abs2, title=:default, ylim_zoom_bool=false, size=(800,600),ylabel="",label=false)
    Omega = operation.(dropdims(sum(result.Omega; dims=3); dims=3))
    z_vec = result.z_vec

    Omega_z = vec(sum(Omega; dims=1))
    Omega_z ./= Omega_z[begin]

    min, max = extrema(Omega_z)
    ymin = (all(>(0), Omega_z)) ? 0.0 : min-0.05abs(min)
    ymin = (ylim_zoom_bool) ? min-0.05abs(min) : ymin
    ymax = max+0.05abs(max)
    ylims = (ymin, ymax)
    x_axis = z_vec
    xlims = extrema(x_axis)

    title = (title==:default) ? "sum $(nameof(operation)) Omega over y and t" : title

    return plot(x_axis, Omega_z; xlabel="\$z\$", ylabel=ylabel, title=title, label=label, ylims=ylims, xlims=xlims, size=size)
end

function plot_2d_superimposed_z_ysum_line(result; operation=abs2, nslices=10, title=:default, zoom=true)
    Omega = operation.(dropdims(sum(result.Omega; dims=3); dims=3))
    Omega ./= maximum(Omega)

    t_vec, z_vec = result.time_vec, result.z_vec
    title = (title==:default) ? "sum $(nameof(operation)) Omega over y" : title
    ylims = extrema(Omega)
    x_axis = t_vec
    xlims = extrema(x_axis)

    if zoom
        startindex,endindex = zoomingin(Omega)
        xlims = (t_vec[startindex], t_vec[endindex])
    end

    colgrad = cgrad([:blue, :yellow, :red])
    colors = colgrad[range(0, 1, length=nslices)]
    fig = plot(; xlabel="\$t\$", title=title, ylims=ylims, xlims=xlims)

    firstlast = i -> (i == 1) ? "first" : ((i==nslices) ? "last" : ((i==nslices÷2+1) ? "middle" : false))

    delta = floor(Int, length(z_vec) / nslices)


    for i in 1:nslices
        plot!(x_axis, Omega[:, i*delta-(delta-1)], c=colors[i], label=firstlast(i))
    end

    return fig
end


function plot_2d_z_yslices_line(result; operation=abs2, nslices=10, title=:default)
    Omega = operation.(result.Omega[:, :, (end÷2+1):end])
    z_vec = result.z_vec
    title = (title==:default) ? "sum $(nameof(operation)) Omega and t" : title
    Omega_z = dropdims(sum(Omega; dims=1); dims=1)
    Omega_z ./= maximum(Omega_z)
    min, max = extrema(Omega_z)
    ylims = (0.0, max+0.05abs(max))
    x_axis = z_vec
    xlims = extrema(x_axis)

    colgrad = cgrad([:blue, :yellow, :red])
    colors = colgrad[range(0, 1, length=nslices)]
    fig = plot(; xlabel="\$z\$", title=title, ylims=ylims, xlims=xlims, legend=false)

    delta = floor(Int, length(result.y_vec) / 2nslices)


    for i in 1:nslices
        plot!(x_axis, Omega_z[:, i*delta-(delta-1)], c=colors[i])
    end

    return fig
end

function plot_2d_zy_pulseprofile_heatmap(result; operation=abs2, title=:default)
    Omega = operation.(result.Omega)
    Omega ./= maximum(Omega)
    y_vec, z_vec = result.y_vec, result.z_vec

    tindex = argmax(Omega[:, end÷2+1, end÷2+1])
    Omega_tslice = Omega[tindex, :, :]
    Omega_tslice ./= maximum(Omega_tslice)

    title = (title==:default) ? "$(nameof(operation)) Omega time slice" : title

    return heatmap(y_vec, z_vec, Omega_tslice; c=:viridis, xlabel="\$y\$", ylabel="\$z\$", title=title, cbar=false)
end


function plot_2d_superimposed_lines_against_t(result; operation=abs2, nslices=10, title=:default, zoom=false)
    Omega = operation.(dropdims(sum(result.Omega; dims=3); dims=3))
    Omega ./= maximum(Omega)

    t_vec, z_vec = result.time_vec, result.z_vec
    title = (title==:default) ? "$(nameof(operation)) Omega" : title
    ylims = extrema(Omega)
    x_axis = t_vec
    xlims = extrema(x_axis) #watch it

    if zoom
        startindex,endindex = zoomingin(Omega)
        xlims = (t_vec[startindex], t_vec[endindex])
    end

    colgrad = cgrad([:blue, :yellow, :red])
    colors = colgrad[range(0, 1, length=nslices)]
    fig = plot(; xlabel="\$t\$", title=title, ylims=ylims, xlims=xlims)

    firstlast = i -> (i == 1) ? "first" : ((i==nslices) ? "last" : ((i==nslices÷2+1) ? "middle" : false))

    delta = floor(Int, length(z_vec) / nslices)

    for i in 1:nslices
        plot!(x_axis, Omega[:, i*delta-(delta-1)], c=colors[i], label=firstlast(i))
    end

    return fig
end

function plot_2sd_echo_efficiency(result)
    Omega = abs2.(dropdims(sum(result.Omega; dims=3); dims=3))

    z_vec = result.z_vec
    t_vec = result.time_vec

    input = result.cfg.pulses[1][1]
    retrieval = result.cfg.pulses[2][1]
    tau = retrieval.center - input.center
    echo_time = input.center + 2tau

    input_peaks = zeros(length(z_vec))
    echo_peaks = zeros(length(z_vec))

    tau_index = searchsortedfirst(t_vec, tau)
    input_index = searchsortedfirst(t_vec, input.center)
    echo_index = searchsortedfirst(t_vec, echo_time)
    width_index = searchsortedfirst(t_vec, input.width)

    inputrange = clamp.((input_index-3width_index):(input_index+3width_index), 1, length(t_vec))

    for i in 1:length(z_vec)
        input_peaks[i] = sum(Omega[inputrange, i])
        echo_peaks[i] = sum(Omega[(echo_index-2width_index):end, i])
    end

    efficiency = echo_peaks ./ input_peaks[1]

    fig = plot(z_vec, efficiency; title="Echo efficiency", xlabel="\$z\$")

    return fig
end

function plot_2d_superimposed_lines_against_y(result; operation=abs2, nslices=10, title=:default, size=(800,600), ylabel="")
    Omega = operation.(result.Omega[end÷2+1,:,:])
    Omega ./= maximum(Omega)

    y_vec = result.y_vec
    z_vec = result.z_vec
    title = (title==:default) ? "$(nameof(operation)) Omega" : title
    ylims = extrema(Omega)
    x_axis = y_vec
    xlims = extrema(x_axis) #watch it

    colgrad = cgrad([:blue, :yellow, :red])
    colors = colgrad[range(0, 1, length=nslices)]
    fig = plot(; xlabel="\$y\$", ylabel=ylabel, title=title, ylims=ylims, xlims=xlims,size=size)

    firstlast = i -> (i == 1) ? "first" : ((i==nslices) ? "last" : ((i==nslices÷2+1) ? "middle" : false))

    indices = round.(Int, range(1, length(z_vec), length=nslices))

    for (i, j) in enumerate(indices)
        plot!(x_axis, Omega[j, :], c=colors[i], label=firstlast(i))
    end

    return fig
end

function plots_for_jevon_2sd_echo(result)
    heatmap = plot_2d_tz_ysum_heatmap(result; operation=abs2, title="Pulse intensity")
    superimposed_zlines_against_t = plot_2d_superimposed_lines_against_t(result; title="Pulse intensity slices in \$z\$", zoom=true, nslices=10)
    echo_efficiency = plot_2sd_echo_efficiency(result)

    return plot(
        heatmap, superimposed_zlines_against_t, echo_efficiency,
        layout=(1, 3), size=(1200, 400), margins=2mm, framestyle=:box, left_margin=4mm, bottom_margin=4mm)
end

function plots_for_jevon_2sd_SIT(result)
    heatmap = plot_2d_tz_ysum_heatmap(result; operation=abs2, title="Pulse intensity")
    #area = plot_2d_z_ysum_line(result; operation=real, title="Total pulse area")
    pulse_profile_heatmap = plot_2d_zy_pulseprofile_heatmap(result; operation=abs2, title="Pulse intensity profile time slice")
    energy = plot_2d_z_ysum_line(result; operation=abs2, title="Total pulse energy")
    superimposed = plot_2d_superimposed_z_ysum_line(result; operation=abs2, nslices=20, title="Pulse intensity slices in z",zoom=true)
    superimposed_area_yslices = plot_2d_z_yslices_line(result; operation=real, nslices=30, title="Pulse area slices in y against z")
    superimposed_intensity_yslices = plot_2d_z_yslices_line(result; operation=abs2, nslices=30, title="Pulse intensity slices in y against z")

    return plot(
        heatmap, pulse_profile_heatmap, superimposed_area_yslices,
        superimposed, energy, superimposed_intensity_yslices,
        layout=(2, 3), size=(1600, 1000), margins=3mm, framestyle=:box, left_margin=5mm, bottom_margin=5mm)
end


#-------------------------------------------

function plot_1d_tz_heatmap(result; operation=abs2, title=:default)
    Omega = operation.(dropdims(result.Omega; dims=3))
    Omega ./= maximum(Omega)
    t_vec, z_vec = result.time_vec, result.z_vec
    title = (title==:default) ? "$(nameof(operation)) Omega" : title
    return heatmap(t_vec, z_vec, transpose(Omega[:, :]);
        c=:viridis, xlabel="\$t\$", ylabel="\$z\$", title=title)
end

function plot_1d_superimposed_lines_against_t(result; operation=abs2, nslices=10, title=:default, zoom=false, size=(800,600), ylabel="")
    Omega = operation.(dropdims(result.Omega; dims=3))
    Omega ./= maximum(Omega)

    t_vec, z_vec = result.time_vec, result.z_vec
    title = (title==:default) ? "$(nameof(operation)) Omega" : title
    ylims = extrema(Omega)
    x_axis = t_vec
    xlims = extrema(x_axis) #watch it

    if zoom
        startindex,endindex = zoomingin(Omega)
        xlims = (t_vec[startindex], t_vec[endindex])
    end

    colgrad = cgrad([:blue, :yellow, :red])
    colors = colgrad[range(0, 1, length=nslices)]
    fig = plot(; xlabel="\$t\$", ylabel=ylabel, title=title, ylims=ylims, xlims=xlims,size=size)

    firstlast = i -> (i == 1) ? "first" : ((i==nslices) ? "last" : ((i==nslices÷2+1) ? "middle" : false))

    indices = round.(Int, range(1, length(z_vec), length=nslices))

    for (i, j) in enumerate(indices)
        plot!(x_axis, Omega[:, j], c=colors[i], label=firstlast(i))
    end

    return fig
end

function plot_1d_superimposed_lines_against_z(result; operation=abs2, nslices=10, title=:default)
    Omega = operation.(dropdims(result.Omega; dims=3))
    Omega ./= maximum(Omega)

    t_vec, z_vec = result.time_vec, result.z_vec
    title = (title==:default) ? "$(nameof(operation)) Omega" : title
    min, max = extrema(Omega)
    ylims = (0.0, max+0.05abs(max))
    x_axis = z_vec
    xlims = extrema(x_axis) #watch it

    colgrad = cgrad([:blue, :yellow, :red])
    colors = colgrad[range(0, 1, length=nslices)]
    fig = plot(; xlabel="z", title=title, ylims=ylims, xlims=xlims)

    firstlast = i -> (i == 1) ? "first" : ((i==nslices) ? "last" : ((i==nslices÷2+1) ? "middle" : false))

    delta = floor(Int, length(t_vec) / nslices)


    for i in 1:nslices
        plot!(x_axis, Omega[i*delta-(delta-1), :], c=colors[i], label=firstlast(i))
    end

    return fig
end

function plot_1sd_echo_efficiency(result)
    Omega = abs2.(dropdims(result.Omega; dims=3))

    z_vec = result.z_vec
    t_vec = result.time_vec

    input = result.cfg.pulses[1][1]
    retrieval = result.cfg.pulses[2][1]
    tau = retrieval.center - input.center
    echo_time = input.center + 2tau

    input_peaks = zeros(length(z_vec))
    echo_peaks = zeros(length(z_vec))

    tau_index = searchsortedfirst(t_vec, tau)
    input_index = searchsortedfirst(t_vec, input.center)
    echo_index = searchsortedfirst(t_vec, echo_time)
    width_index = searchsortedfirst(t_vec, input.width)

    inputrange = clamp.((input_index-3width_index):(input_index+3width_index), 1, length(t_vec))

    for i in 1:length(z_vec)
        input_peaks[i] = sum(Omega[inputrange, i])
        echo_peaks[i] = sum(Omega[(echo_index-2width_index):end, i])
    end

    efficiency = echo_peaks ./ input_peaks[1]

    fig = plot(z_vec, efficiency; title="Echo efficiency", xlabel="\$z\$")

    return fig
end

function plots_for_jevon_1sd_echo(result)
    heatmap = plot_1d_tz_heatmap(result; operation=abs2, title="Pulse intensity")
    superimposed_zlines_against_t = plot_1d_superimposed_lines_against_t(result; title="Pulse intensity slices in \$z\$", zoom=true, nslices=10)
    echo_efficiency = plot_1sd_echo_efficiency(result)

    return plot(
        heatmap, superimposed_zlines_against_t, echo_efficiency,
        layout=(1, 3), size=(1600, 600), margins=2mm, framestyle=:box, left_margin=4mm, bottom_margin=4mm)
end

function plots_for_jevon_1sd_SIT(result)
    heatmap = plot_1d_tz_heatmap(result; operation=abs2, title="Pulse intensity")
    superimposed_zlines_against_t = plot_1d_superimposed_lines_against_t(result; title="Pulse intensity slices in \$z\$", zoom=true, nslices=10)
    energy = plot_2d_z_ysum_line(result; operation=abs2, title="Sum of pulse intensity across \$t\$", ylim_zoom_bool=true)
    area = plot_2d_z_ysum_line(result; operation=real, title="Pulse area", ylim_zoom_bool=true)
    #superimposed_tlines_against_z = plot_1d_superimposed_lines_against_z(result; operation=abs2, nslices=10, title="Pulse intensity slices in \$t\$")

    return plot(
        heatmap, energy,
        superimposed_zlines_against_t, area,
        layout=(2, 2), size=(800, 500), margins=2mm, framestyle=:box, left_margin=4mm, bottom_margin=4mm)
end

#------------------------

function plot_0sd_polarisation(result)
    P = abs2.(result.P)

    t_vec = result.time_vec
    #title = (title==:default) ? "" : title
    fig = plot(; ylabel="\$|\\mathcal{P}|^2\$", title="Squared magnitude of the polarisation density against time")

    plot!(t_vec, P; label=false)

    return fig
end

function plot_0sd_pol_complex_parts(result)
    Pu = 2*real(result.P)
    Pv = -2*imag(result.P)

    t_vec = result.time_vec

    fig = plot(; title="Ensemble averaged Bloch components against time")

    plot!(t_vec, Pu; label="averaged u")
    plot!(t_vec, Pv; label="averaged v")

    return fig
end

function plots_for_jevon_0sd(result)
    pol_abs2 = plot_0sd_polarisation(result)
    pol_real_n_imag = plot_0sd_pol_complex_parts(result)


    return plot(
        pol_abs2, pol_real_n_imag; layout=(2, 1), size=(800, 600)
    )
end

#--------------------------------------------

using LsqFit

function plot_fitted_SIT(result;zoom=true)
    Omega = abs2.(dropdims(sum(result.Omega;dims=3); dims=3))
    Omega ./= maximum(Omega)
    t_vec = result.time_vec

    sech2fit(x, p)  = @. p[1] / cosh((x-p[2])/p[3])^2 + p[4]
    gaussianfit(x, p) = @. p[1] * exp(-0.5*((x-p[2])/p[3])^2) + p[4]

    # Your simulation data
    x =  t_vec
    y1 = Omega[:,begin]
    y2 = Omega[:,end]

    # Initial guesses
    p01 = [maximum(y1), x[argmax(y1)], (maximum(x)-minimum(x))/10, 0.0]
    p02 = [maximum(y2), x[argmax(y2)], (maximum(x)-minimum(x))/10, 0.0]

    # Fit both
    fit_g = curve_fit(gaussianfit, x, y1, p01)
    fit_s = curve_fit(sech2fit, x, y2, p02)

    # Best-fitting parameters
    println(fit_g.param)
    println(fit_s.param)

    fig = plot_1d_superimposed_lines_against_t(result; title="", zoom=true, nslices=10,size=(800,600),ylabel="\$I(z)/I_{max}\$")
    plot!(x, gaussianfit(x, fit_g.param), label="Gaussian fit", lc=:black, lw=1.5,ls=:dash)
    plot!(x, sech2fit(x, fit_s.param), label="Sech^2 fit", lc=:black, lw=1.5,ls=:dot)

    if zoom
        startindex,endindex = zoomingin(Omega)
        xlims = (t_vec[startindex], t_vec[endindex])
    end
    difffig = plot(x, y1 - gaussianfit(x, fit_g.param); label="Input face data minus Gaussian fit", size=(800,600),xlims=xlims)
    plot!(x, y2 - sech2fit(x, fit_s.param); label="Output face data minus sech^2 fit")
    

    fig = plot(fig, difffig; layout=(1,2), size=(1600,600), bottom_margin=10mm, left_margin=10mm)

    return fig
end

function plot_fitted_attenuation(result)
    Omega = abs2.(dropdims(sum(result.Omega;dims=3); dims=3))
    Omega ./= maximum(Omega)
    t_vec = result.time_vec
    z_vec = result.z_vec

    # Your simulation data
    x = z_vec
    y = vec(sum(Omega[:,:],dims=1))
    y ./= maximum(y)

    z0 = x[begin]
    exponentialfit(x, p) = @. p[1] * exp(-(x-z0)/p[2]) + p[3]

    # Initial guesses
    p0 = [maximum(y), (maximum(x)-minimum(x))/10, 0.0]

    # Fit both
    fit_e = curve_fit(exponentialfit, x, y, p0)
    expfit = exponentialfit(x, fit_e.param)


    fig = plot_2d_z_ysum_line(result; operation=abs2, title="", ylim_zoom_bool=true, size=(800,600), ylabel="\$E(z)/E_{max}\$", label="Data")
    plot!(x, expfit, label="Exponential fit", lc=:black, lw=1.5,ls=:dash)

    difffig = plot(x, y - expfit; label="Data minus exponential fit", size=(800,600), xlabel="\$z\$", ylabel="\$E(z)/E_{max}\$")
    
    fig = plot(fig, difffig; layout=(1,2), size=(800,300), bottom_margin=2mm, left_margin=3mm)

    return fig
end


function plot_fitted_diffraction(result)
    Omega = abs2.(result.Omega[end÷2+1,:,:]) #halfway in time
    Omega ./= maximum(Omega)
    y_vec = result.y_vec

    gaussianfit(x, p) = @. p[1] * exp(-0.5*((x-p[2])/p[3])^2) + p[4]

    # Your simulation data
    x = y_vec
    y1 = Omega[begin,:]
    y2 = Omega[end,:]

    # Initial guesses
    p01 = [maximum(y1), x[argmax(y1)], (maximum(x)-minimum(x))/10, 0.0]
    p02 = [maximum(y2), x[argmax(y2)], (maximum(x)-minimum(x))/10, 0.0]

    # Fit both
    fit_g1 = curve_fit(gaussianfit, x, y1, p01)
    fit_g2 = curve_fit(gaussianfit, x, y2, p02)

    # Best-fitting parameters
    println(fit_g1.param)
    println(fit_g2.param)

    fig = plot_2d_superimposed_lines_against_y(result; title="", nslices=10,size=(800,600),ylabel="\$I(z)/I_{max}\$")
    plot!(x, gaussianfit(x, fit_g1.param), label="Fit for input face", lc=:black, lw=1.5,ls=:dash)
    plot!(x, gaussianfit(x, fit_g2.param), label="Fit for output face", lc=:black, lw=1.5,ls=:dot)

    difffig = plot(x, y1 - gaussianfit(x, fit_g1.param); label="Input face data minus fit", size=(800,600))
    plot!(x, y2 - gaussianfit(x, fit_g2.param); label="Output face data minus fit")
    

    fig = plot(fig, difffig; layout=(1,2), size=(1600,600), bottom_margin=10mm, left_margin=10mm)

    return fig
end