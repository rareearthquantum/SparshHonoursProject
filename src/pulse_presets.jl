# Gaussian pulse sequences for the example experiments.
# Centres and temporal widths are fractions of Tf - Ti.
# These functions construct pulses only; they do not select propagation parameters.

function default_echo_2d_pulses(;Ti::Real=0.0, Tf::Real=1.0, y_pulse_width::Real=1.0)
    duration = Tf - Ti
    return [
        (PulseParams(Ti + duration/10, duration/25, pi/2), PulseParams(0.0, y_pulse_width, 1.0)),
        (PulseParams(Ti + 4duration/10, duration/100, pi), PulseParams(0.0, y_pulse_width, 1.0))
    ]
end

function default_soliton_2d_pulses(;Ti::Real=0.0, Tf::Real=1.0, y_pulse_width::Real=1.0)
    duration = Tf - Ti
    return [
        (PulseParams(Ti + 2duration/10, duration/50, 2pi), PulseParams(0.0, y_pulse_width, 1.0))
    ]
end

function default_4pi_2d_pulses(;Ti::Real=0.0, Tf::Real=1.0, y_pulse_width::Real=1.0)
    duration = Tf - Ti
    return [
        (PulseParams(Ti + 2duration/10, duration/20, 4pi), PulseParams(0.0, y_pulse_width, 1.0))
    ]
end

function default_6pi_2d_pulses(;Ti::Real=0.0, Tf::Real=1.0, y_pulse_width::Real=1.0)
    duration = Tf - Ti
    return [
        (PulseParams(Ti + 2duration/10, duration/20, 6pi), PulseParams(0.0, y_pulse_width, 1.0))
    ]
end

function default_weak_2d_pulses(;Ti::Real=0.0, Tf::Real=1.0, y_pulse_width::Real=1.0)
    duration = Tf - Ti
    return [
        (PulseParams(Ti + duration/2, duration/10, pi/8), PulseParams(0.0, y_pulse_width, 1.0))
    ]
end

function default_pi_2d_pulses(;Ti::Real=0.0, Tf::Real=1.0, y_pulse_width::Real=1.0)
    duration = Tf - Ti
    return [
        (PulseParams(Ti + 3duration/10, duration/50, pi), PulseParams(0.0, y_pulse_width, 1.0))
    ]
end
