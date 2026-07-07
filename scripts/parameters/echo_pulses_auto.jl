pulses = PulseParams[
    PulseParams(center=1u"μs", width=0.2u"μs", area=pi/8),
    PulseParams(center=4u"μs", width=0.1u"μs", area=pi),
]

# Ti, Tf, Nt, d_width, and Nd are calculated from the pulses.
# Zf, Nz, alpha, samples_per_width, and padding_widths can be overridden here.
config_from_pulses(pulses; Zf=10u"cm", alpha=100u"s^-1*m^-1")
