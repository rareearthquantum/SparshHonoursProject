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
