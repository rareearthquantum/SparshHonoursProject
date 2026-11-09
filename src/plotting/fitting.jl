using LsqFit

# Curve-fit comparisons for the existing research plots.
# Shared plotting helpers are defined in helpers.jl.

"""
    _gaussian_fit(x, p)

Evaluate a Gaussian plus constant offset at `x`, with
`p = [amplitude, center, standard_deviation, offset]`.
"""
_gaussian_fit(x, p) = @. p[1] * exp(-0.5 * ((x-p[2])/p[3])^2) + p[4]
"""
    _sech_squared_fit(x, p)

Evaluate a squared hyperbolic secant plus constant offset at `x`, with
`p = [amplitude, center, width, offset]`.
"""
_sech_squared_fit(x, p) = @. p[1] / cosh((x-p[2])/p[3])^2 + p[4]

"""
    _fit_initial_parameters(x, y)

Return initial fit parameters `[maximum(y), x[argmax(y)], span(x)/10, 0]`
for Gaussian or squared-sech fitting. `x` and `y` must be nonempty
corresponding sample vectors.
"""
function _fit_initial_parameters(x, y)
    return [maximum(y), x[argmax(y)], (maximum(x)-minimum(x))/10, 0.0]
end

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

"""
    plot_fitted_attenuation(result)

Sum intensity over time and transverse position, normalise by its
maximum, and fit an exponential decay plus offset along z. Return a
two-panel plot showing the data and fit alongside their residuals.
"""
function plot_fitted_attenuation(result)
    field_values = _field_2d_summed_over_y(result, abs2)
    _normalise_by_max!(field_values)
    z_vec = result.z_vec

    # Your simulation data
    x = z_vec
    y = vec(sum(field_values[:, :], dims=1))
    _normalise_by_max!(y)

    z0 = x[begin]
    # Exponential model: p = [amplitude, decay length, offset], anchored at z0.
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
