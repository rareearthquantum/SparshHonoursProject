"""
Evaluate an area-normalised Gaussian at `u`. `center` is the peak position,
`width` is the positive standard deviation, and `area` is its integral over
the real line.
"""
function pulse(u, center, width, area)
    return area / (sqrt(2pi) * width) * exp(-0.5 * ((u - center) / width)^2)
end

"""
Evaluate a unit-peak Gaussian with the given centre and positive standard
deviation. The four-argument method accepts `area` for compatibility but
ignores it.
"""
function gaussian(u, center, width)
    return exp(-0.5 * ((u - center) / width)^2)
end

"""
Evaluate a unit-peak Gaussian with the given centre and positive standard
deviation. The four-argument method accepts `area` for compatibility but
ignores it.
"""
gaussian(u, center, width, area) = exp(-0.5 * ((u - center) / width)^2)

"""
Evaluate the area-normalised Gaussian at `t` using `params.center`,
`params.width` (positive standard deviation), and `params.area`.
"""
pulse(t, params::PulseParams) = pulse(t, params.center, params.width, params.area)

"""
Evaluate a unit-peak Gaussian at `u` using `pulseparams.center` and
`pulseparams.width` as the positive standard deviation. The area field
is ignored.
"""
gaussian(u, pulseparams::PulseParams) = gaussian(u, pulseparams.center, pulseparams.width, pulseparams.area)
