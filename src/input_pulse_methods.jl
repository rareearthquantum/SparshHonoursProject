function pulse(
    u,
    center,
    width,
    area
)
    return area * (1 / (sqrt(2pi) * width)) * exp(-0.5 * ((u - center) / width)^2)
end

function pulse(
    u,
    pulse_params
)
    center, width, area = pulse_params
    pulse(u, center, width, area)
end


struct PulseParams
    center
    width
    area
end

function initialise_pulses(t_grid, y_grid, pulse_params)
    tpp = (0.0, 0.0, 0.0)
    ypp = (0.0, 0.0, 0.0)
    pulses = zeros(Complex, length(t_grid), length(y_grid))
    for i in eachindex(pulse_params)
        tpp = pulse_params[i][1]
        ypp = pulse_params[i][2]
        pulses .+= pulse.(t_grid, tpp.center, tpp.width, tpp.area) .* pulse.(y_grid, ypp.center, ypp.width, ypp.area)'
    end
    return pulses
end


function flat_pulse(u,ui,uf)::Real
    return (ui <= u && u <= uf)
end