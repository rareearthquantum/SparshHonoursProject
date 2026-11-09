"""
    plot_fitted_diffraction(result)

At the middle time sample, fit Gaussian transverse intensity profiles
at the input and output planes after shared maximum normalisation.
Print fitted parameters and return a two-panel plot of profiles, fits,
and residuals.
"""
function plot_fitted_diffraction(result)
    field_values = abs2.(result.Omega[end÷2+1, :, :]) #halfway in time
    _normalise_by_max!(field_values)
    y_vec = result.y_vec

    # Your simulation data
    x = y_vec
    y1 = field_values[begin, :]
    y2 = field_values[end, :]

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
