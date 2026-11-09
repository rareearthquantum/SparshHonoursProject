# Gaussian pulse sequences for the example experiments.
# Centres and temporal widths are fractions of Tf - Ti.
# These functions construct pulses only; they do not select propagation parameters.

"""
    default_echo_2d_pulses(;Ti::Real=0.0, Tf::Real=1.0, y_pulse_width::Real=1.0)

Return a vector of `(temporal, transverse)` `PulseParams` pairs for
a pi/2 input pulse at 10% of the time window (width 4%) and a pi
retrieval pulse at 40% (width 1%). Each transverse profile is a unit-peak Gaussian
centred at zero with standard deviation `y_pulse_width`. Requires
`Tf > Ti` and positive transverse width; propagation parameters are
chosen separately.
"""
function default_echo_2d_pulses(;Ti::Real=0.0, Tf::Real=1.0, y_pulse_width::Real=1.0)
    duration = Tf - Ti
    return [
        (PulseParams(Ti + duration/10, duration/25, pi/2), PulseParams(0.0, y_pulse_width, 1.0)),
        (PulseParams(Ti + 4duration/10, duration/100, pi), PulseParams(0.0, y_pulse_width, 1.0))
    ]
end

"""
    default_soliton_2d_pulses(;Ti::Real=0.0, Tf::Real=1.0, y_pulse_width::Real=1.0)

Return a vector of `(temporal, transverse)` `PulseParams` pairs for
a single 2pi temporal pulse at 20% of the time window
with standard deviation 2% of its duration. Each transverse profile is a unit-peak Gaussian
centred at zero with standard deviation `y_pulse_width`. Requires
`Tf > Ti` and positive transverse width; propagation parameters are
chosen separately.
"""
function default_soliton_2d_pulses(;Ti::Real=0.0, Tf::Real=1.0, y_pulse_width::Real=1.0)
    duration = Tf - Ti
    return [
        (PulseParams(Ti + 2duration/10, duration/50, 2pi), PulseParams(0.0, y_pulse_width, 1.0))
    ]
end

"""
    default_4pi_2d_pulses(;Ti::Real=0.0, Tf::Real=1.0, y_pulse_width::Real=1.0)

Return a vector of `(temporal, transverse)` `PulseParams` pairs for
a single 4pi temporal pulse at 20% of the time window
with standard deviation 5% of its duration. Each transverse profile is a unit-peak Gaussian
centred at zero with standard deviation `y_pulse_width`. Requires
`Tf > Ti` and positive transverse width; propagation parameters are
chosen separately.
"""
function default_4pi_2d_pulses(;Ti::Real=0.0, Tf::Real=1.0, y_pulse_width::Real=1.0)
    duration = Tf - Ti
    return [
        (PulseParams(Ti + 2duration/10, duration/20, 4pi), PulseParams(0.0, y_pulse_width, 1.0))
    ]
end

"""
    default_6pi_2d_pulses(;Ti::Real=0.0, Tf::Real=1.0, y_pulse_width::Real=1.0)

Return a vector of `(temporal, transverse)` `PulseParams` pairs for
a single 6pi temporal pulse at 20% of the time window
with standard deviation 5% of its duration. Each transverse profile is a unit-peak Gaussian
centred at zero with standard deviation `y_pulse_width`. Requires
`Tf > Ti` and positive transverse width; propagation parameters are
chosen separately.
"""
function default_6pi_2d_pulses(;Ti::Real=0.0, Tf::Real=1.0, y_pulse_width::Real=1.0)
    duration = Tf - Ti
    return [
        (PulseParams(Ti + 2duration/10, duration/20, 6pi), PulseParams(0.0, y_pulse_width, 1.0))
    ]
end

"""
    default_weak_2d_pulses(;Ti::Real=0.0, Tf::Real=1.0, y_pulse_width::Real=1.0)

Return a vector of `(temporal, transverse)` `PulseParams` pairs for
a single pi/8 temporal pulse at the midpoint of the time window
with standard deviation 10% of its duration. Each transverse profile is a unit-peak Gaussian
centred at zero with standard deviation `y_pulse_width`. Requires
`Tf > Ti` and positive transverse width; propagation parameters are
chosen separately.
"""
function default_weak_2d_pulses(;Ti::Real=0.0, Tf::Real=1.0, y_pulse_width::Real=1.0)
    duration = Tf - Ti
    return [
        (PulseParams(Ti + duration/2, duration/10, pi/8), PulseParams(0.0, y_pulse_width, 1.0))
    ]
end

"""
    default_pi_2d_pulses(;Ti::Real=0.0, Tf::Real=1.0, y_pulse_width::Real=1.0)

Return a vector of `(temporal, transverse)` `PulseParams` pairs for
a single pi temporal pulse at 30% of the time window
with standard deviation 2% of its duration. Each transverse profile is a unit-peak Gaussian
centred at zero with standard deviation `y_pulse_width`. Requires
`Tf > Ti` and positive transverse width; propagation parameters are
chosen separately.
"""
function default_pi_2d_pulses(;Ti::Real=0.0, Tf::Real=1.0, y_pulse_width::Real=1.0)
    duration = Tf - Ti
    return [
        (PulseParams(Ti + 3duration/10, duration/50, pi), PulseParams(0.0, y_pulse_width, 1.0))
    ]
end
