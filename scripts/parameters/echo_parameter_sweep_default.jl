pulses = PulseParams[
    PulseParams(center=1u"μs", width=0.4u"μs", area=pi/8),
    PulseParams(center=4u"μs", width=0.1u"μs", area=pi)
]

(
    base_config=config_from_pulses(pulses; Zf=10u"m", Nz=256),
    sweep_parameters=(
        alpha=[1.0, 5.0, 10.0, 20.0]u"s^-1*m^-1",
    ),
)
