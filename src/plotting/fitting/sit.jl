"""
    plot_fitted_SIT(result; zoom=true)

Fit the transversely summed input intensity to a Gaussian and output
intensity to a squared sech, after shared maximum normalisation. Print
the fitted parameters and return a two-panel plot of data, fits, and
residuals. `zoom=true` crops time to the above-threshold field region.
"""
function plot_fitted_SIT(result; zoom=true)
    field_values = _field_2d_summed_over_y(result, abs2)
    _normalise_by_max!(field_values)
    t_vec = result.t_vec

    # Your simulation data
    x = t_vec
    y1 = field_values[:, begin]
    y2 = field_values[:, end]

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
    xlims = _plot_xlims(t_vec, field_values, zoom)

    fig = plot(x, y1; title="", ylabel="\$I/I_{max}\$", xlabel="\$t\$", label="Input face data", xlims=xlims, legend=:topleft, st=:scatter, shape=:circle, ms=3, mc=:red, msc=:red)
    plot!(x, y2; label="Output face data", st=:scatter, shape=:circle, ms=3, mc=:blue, msc=:blue)
    plot!(x, gaussian_curve, label="Gaussian fit", st=:line, lc=:red, ls=:solid)
    plot!(x, sech_curve, label="Sech^2 fit", st=:line, lc=:blue, ls=:solid)

    difffig = plot(x, y1 - gaussian_curve; label="Input face data-fit difference", xlims=xlims, ylabel="\$I/I_{max}\$", xlabel="\$t\$", legend=:topleft)
    plot!(x, y2 - sech_curve; label="Output face data-fit difference")

    return _compose_plots(fig, difffig; layout=(1, 2))
end
