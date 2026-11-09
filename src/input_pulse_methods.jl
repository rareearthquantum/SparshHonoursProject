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

struct PulseParams
    center::Float64
    width::Float64
    area::Float64
end

"""
Construct Gaussian or top-hat pulse parameters, converting keyword values
to `Float64`. For Gaussians, `width` is the standard deviation; for top hats,
it is the full support width. Validation is performed by `EchoConfig`.
"""
PulseParams(; center, width, area) =
    PulseParams(Float64(center), Float64(width), Float64(area))

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

"""
Return the sum of area-normalised temporal Gaussians at time `t`, using
the supplied vector of `PulseParams`.
"""
pulse_sum(t, pulses::AbstractVector{<:PulseParams}) = sum(params -> pulse(t, params), pulses)

"""
Return a `ComplexF64` matrix indexed by `(time, transverse position)`.
Each entry of `pulse_params` is a pair of splattable `(center, width, area)`
tuples; both temporal and transverse Gaussians are area-normalised here.
"""
function initialise_pulses(t_grid, y_grid, pulse_params)
    pulses = zeros(ComplexF64, length(t_grid), length(y_grid))
    for (t_params, y_params) in pulse_params
        pulses .+= pulse.(t_grid, t_params...) .* pulse.(y_grid, y_params...)'
    end
    return pulses
end

"""
Evaluate a separable pulse at `(t, y)` from a `(temporal, transverse)` pair
of `PulseParams`. The temporal Gaussian has the specified area; the
transverse Gaussian has unit peak and ignores its area field.
"""
function pulse_2d(t, y, ty_param)
    return pulse(t,ty_param[1])*gaussian(y,ty_param[2])
end

"""
Return the summed Gaussian field at `(t, y)` for a vector of
`(temporal, transverse)` `PulseParams` pairs. See [`pulse_2d`](@ref).
"""
pulse_2d_sum(t, y, pulses_2d::AbstractVector{NTuple{2,PulseParams}}) = sum(params -> pulse_2d(t, y, params), pulses_2d)


"""
Return an area-normalised rectangular pulse as `Float64`. The value is
`area / width` when `abs(u - center) <= width / 2`, including the endpoints,
and zero otherwise. `width` must be positive.
"""
top_hat_pulse(u, center, width, area)::Float64 = ( abs(u-center) <= width/2 ) ? area/width : 0
"""
Evaluate an area-normalised top hat at `u` using the centre, positive full
width, and area stored in `pulseparams`. Return `Float64`; the support
includes both endpoints.
"""
top_hat_pulse(u,  pulseparams::PulseParams)::Float64 = top_hat_pulse(u, pulseparams.center, pulseparams.width, pulseparams.area)

"""
Evaluate a temporal top hat times a unit-peak transverse Gaussian at
`(t, y)`, using a `(temporal, transverse)` pair of `PulseParams`.
"""
tophat_2dish_pulse(t, y, ty_param) = top_hat_pulse(t,ty_param[1])*gaussian(y,ty_param[2])
"""
Sum temporal top hats with Gaussian transverse profiles at `(t, y)`
for the supplied `(temporal, transverse)` `PulseParams` pairs.
"""
tophat_pulse_2dish_sum(t, y, pulses_2d::AbstractVector{NTuple{2,PulseParams}}) = sum(params -> tophat_2dish_pulse(t, y, params), pulses_2d)

"""
Evaluate the product of area-normalised temporal and transverse top hats
at `(t, y)`, using a pair of `PulseParams`.
"""
tophat_2d_pulse(t, y, ty_param) = top_hat_pulse(t,ty_param[1])*top_hat_pulse(y,ty_param[2])
"""
Sum separable top-hat pulses at `(t, y)` for the supplied
`(temporal, transverse)` `PulseParams` pairs.
"""
tophat_pulse_2d_sum(t, y, pulses_2d::AbstractVector{NTuple{2,PulseParams}}) = sum(params -> tophat_2d_pulse(t, y, params), pulses_2d)
