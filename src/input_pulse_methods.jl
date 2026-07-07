using Unitful

seconds(x::Real) = Float64(x)
seconds(x::Unitful.AbstractQuantity) = ustrip(Float64, u"s", x)

function pulse(u, center, width, area)
    return area / (sqrt(2pi) * width) * exp(-0.5 * ((u - center) / width)^2)
end

struct PulseParams
    center::Float64
    width::Float64
    area::Float64
end

PulseParams(; center, width, area) =
    PulseParams(seconds(center), seconds(width), Float64(area))

pulse(t, params::PulseParams) = pulse(t, params.center, params.width, params.area)

pulse_sum(t, pulses::AbstractVector{<:PulseParams}) = sum(params -> pulse(t, params), pulses)

function initialise_pulses(t_grid, y_grid, pulse_params)
    pulses = zeros(ComplexF64, length(t_grid), length(y_grid))
    for (t_params, y_params) in pulse_params
        pulses .+= pulse.(t_grid, Ref(t_params)) .* pulse.(y_grid, Ref(y_params))'
    end
    return pulses
end
