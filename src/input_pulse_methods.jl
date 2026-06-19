function pulse(
    u,
    center,
    width,
    area
)
    normalisation_factor = 1 / (sqrt(2pi) * width)
    return area * normalisation_factor * exp(-0.5 * ((u - center) / width)^2)
end

function pulse(
    u,
    pulse_params
)
    center, width, area = pulse_params
    pulse(u, center, width, area)
end


#Preset pulse in t
function Ein_t(t)
    p1c, p1w, p1A = 5.0, 1.0, pi / 1
    p1 = pulse(t, p1c, p1w, p1A)
    #p2 = pulse(t,18.0,1.0,1.0)
    return p1
end

#Preset pulse in y
function Ein_y(y)
    p_c, p_w, p_A = 0.0, 1.0, 1.0
    p = pulse(y, p_c, p_w, p_A)
    return p
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