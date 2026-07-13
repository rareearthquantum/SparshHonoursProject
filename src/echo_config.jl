const ND_TO_NT_RATIO = 0.2

detuning_width(Nt, Ti, Tf) = (Nt - 1) / (Tf - Ti) # 1/dt
detuning_count(Nt) = ceil(Int, ND_TO_NT_RATIO * Nt)
meters(x::Real) = Float64(x)
meters(x::Unitful.AbstractQuantity) = ustrip(Float64, u"m", x)
coupling_si(x::Real) = Float64(x)
coupling_si(x::Unitful.AbstractQuantity) = ustrip(Float64, u"s^-1*m^-1", x)

function default_echo_pulses(Ti::Real, Tf::Real)
    duration = Tf - Ti
    return PulseParams[
        PulseParams(Ti + duration/10, duration/25, pi/2),
        PulseParams(Ti + 4duration/10, duration/100, pi),
    ]
end

Base.@kwdef struct EchoConfig
    Nt::Int = 1000
    Ti::Float64 = 0.0
    Tf::Float64 = 10.0

    d_width::Float64 = detuning_width(max(1000,Nt), Ti, Tf)
    Nd::Int = detuning_count(max(1000,Nt))

    pulses::Vector{PulseParams} = default_echo_pulses(Ti, Tf)

    Nz::Int = 200
    Zi::Float64 = 0.0
    Zf::Float64 = 10.0

    alpha::Float64 = 1.0
end

"""Build a simulation window and grids from an ordered list of echo pulses."""
function config_from_pulses(pulses;
        samples_per_width=10,
        padding_widths=6.0,
        Zf=10.0,
        Nz::Int=200,
        alpha=1.0e2)
    first_pulse = first(pulses)
    last_pulse = last(pulses)
    Ti = min(0.0, minimum(pulse.center - padding_widths*pulse.width for pulse in pulses))
    pulse_end = maximum(pulse.center + padding_widths*pulse.width for pulse in pulses)

    echo_center = 2last_pulse.center - first_pulse.center
    echo_end = echo_center + padding_widths*first_pulse.width
    Tf = max(pulse_end, echo_end)

    dt = minimum(pulse.width for pulse in pulses) / samples_per_width
    Nt = ceil(Int, (Tf - Ti) / dt) + 1

    return EchoConfig(
        Nt=Nt, Ti=Ti, Tf=Tf,
        pulses=collect(pulses),
        Nz=Nz, Zf=meters(Zf), alpha=coupling_si(alpha),
    )
end

make_detunings(cfg::EchoConfig) = LinRange(-cfg.d_width/2, cfg.d_width/2, cfg.Nd)
make_time_grid(cfg::EchoConfig) = LinRange(cfg.Ti, cfg.Tf, cfg.Nt)
make_z_grid(cfg::EchoConfig) = LinRange(cfg.Zi, cfg.Zf, cfg.Nz)
make_omega_input(cfg::EchoConfig) = t -> pulse_sum(t, cfg.pulses)