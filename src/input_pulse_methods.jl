function pulse(
    u::Float64,
    center::Float64,
    width::Float64,
    area::Float64
)
    normalisation_factor = 1 / (sqrt(pi) * width)
    area * normalisation_factor * exp(-((u - center) / width)^2)
end


#Preset pulse in t
function Ein_t(t::Float64)
    p1c, p1w, p1A = 5.0, 1.0, pi / 1
    p1 = pulse(t, p1c, p1w, p1A)
    #p2 = pulse(t,18.0,1.0,1.0)
    return p1
end

#Preset pulse in y
function Ein_y(y::Float64)
    p_c, p_w, p_A = 0.0, 1.0, 1.0
    p = pulse(y, p_c, p_w, p_A)
    return p
end

struct PulseParams
    center::Float64
    width::Float64
    area::Float64
end

function initialise_pulses(t_grid, y_grid, pulse_params::Array{NTuple{2,PulseParams}})
    tpp = (0.0, 0.0, 0.0)
    ypp = (0.0, 0.0, 0.0)
    pulses = zeros(ComplexF64, length(t_grid), length(y_grid))
    for i in eachindex(pulse_params)
        tpp = pulse_params[i][1]
        ypp = pulse_params[i][2]
        pulses .+= pulse.(t_grid, tpp.center, tpp.width, tpp.area) .* fftshift(fft(pulse.(y_grid, ypp.center, ypp.width, ypp.area)))'
    end
    return pulses
end