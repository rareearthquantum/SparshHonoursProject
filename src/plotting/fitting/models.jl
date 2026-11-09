using LsqFit

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
