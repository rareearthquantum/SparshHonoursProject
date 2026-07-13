using Unitful

seconds(x::Real) = Float64(x)
seconds(x::Unitful.AbstractQuantity) = ustrip(Float64, u"s", x)

function pulse(u, center, width, area)
    return area / (sqrt(2pi) * width) * exp(-0.5 * ((u - center) / width)^2)
end

function gaussian(u, center, width)
    return exp(-0.5 * ((u - center) / width)^2)
end
gaussian(u, center, width, area) = exp(-0.5 * ((u - center) / width)^2)

struct PulseParams
    center::Float64
    width::Float64
    area::Float64
end

PulseParams(; center, width, area) =
    PulseParams(seconds(center), seconds(width), Float64(area))

pulse(t, params::PulseParams) = pulse(t, params.center, params.width, params.area)
gaussian(u, pulseparams::PulseParams) = gaussian(u, pulseparams.center, pulseparams.width, pulseparams.area)

pulse_sum(t, pulses::AbstractVector{<:PulseParams}) = sum(params -> pulse(t, params), pulses)

function initialise_pulses(t_grid, y_grid, pulse_params)
    pulses = zeros(ComplexF64, length(t_grid), length(y_grid))
    for (t_params, y_params) in pulse_params
        pulses .+= pulse.(t_grid, t_params...) .* pulse.(y_grid, y_params...)'
    end
    return pulses
end

function pulse_2d(t, y, ty_param)
    return pulse(t,ty_param[1])*gaussian(y,ty_param[2])
end

pulse_2d_sum(t, y, pulses_2d::AbstractVector{NTuple{2,PulseParams}}) = sum(params -> pulse_2d(t, y, params), pulses_2d)


top_hat_pulse(u, center, width, area)::Float64 = ( abs(u-center) <= width/2 ) ? area/width : 0
top_hat_pulse(u,  pulseparams::PulseParams)::Float64 = top_hat_pulse(u, pulseparams.center, pulseparams.width, pulseparams.area)