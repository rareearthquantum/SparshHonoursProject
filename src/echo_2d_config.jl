const ND_TO_NT_RATIO = 0.2

detuning_width(Nt, Ti, Tf) = (Nt - 1) / (Tf - Ti) # 1/dt
detuning_count(Nt) = ceil(Int, ND_TO_NT_RATIO * Nt)
y_width(y_pulse_width, beta, Z_length) = 6*y_pulse_width*sqrt(1+(2*beta*Z_length/y_pulse_width^2)^2)
meters(x::Real) = Float64(x)
meters(x::Unitful.AbstractQuantity) = ustrip(Float64, u"m", x)
coupling_si(x::Real) = Float64(x)
coupling_si(x::Unitful.AbstractQuantity) = ustrip(Float64, u"s^-1*m^-1", x)

function default_echo_2d_pulses(Ti::Real, Tf::Real; y_pulse_width::Real=1.0)
    duration = Tf - Ti
    return [
        (PulseParams(Ti + duration/10, duration/25, pi/2), PulseParams(0.0, y_pulse_width, 1.0)),
        (PulseParams(Ti + 4duration/10, duration/100, pi), PulseParams(0.0, y_pulse_width, 1.0))
    ]
end

function default_soliton_2d_pulses(Ti::Real, Tf::Real; y_pulse_width::Real=1.0)
    duration = Tf - Ti
    return [
        (PulseParams(Ti + 3duration/10, duration/10, 2pi), PulseParams(0.0, y_pulse_width, 1.0))
    ]
end

Base.@kwdef struct EchoConfig
    Nt::Int = 256
    Ti::Float64 = 0.0
    Tf::Float64 = 10.0

    d_width::Float64 = detuning_width(Nt, Ti, Tf)
    Nd::Int = detuning_count(Nt)

    Nz::Int = 64
    Zi::Float64 = 0.0
    Zf::Float64 = 10.0

    alpha::Float64 = 1.0e1
    beta::Float64 = 1e-2

    Ny::Int = 64
    y_pulse_width::Float64 = 1.0
    Yi::Float64 = -y_width(y_pulse_width, beta, Zf-Zi)/2
    Yf::Float64 = y_width(y_pulse_width, beta, Zf-Zi)/2

    pulses::Vector{NTuple{2,PulseParams}} = default_echo_2d_pulses(Ti, Tf; y_pulse_width)
end

make_detunings(cfg::EchoConfig) = LinRange(-cfg.d_width/2, cfg.d_width/2, cfg.Nd)
make_time_grid(cfg::EchoConfig) = LinRange(cfg.Ti, cfg.Tf, cfg.Nt)
make_z_grid(cfg::EchoConfig) = LinRange(cfg.Zi, cfg.Zf, cfg.Nz)
make_y_grid(cfg::EchoConfig) = LinRange(cfg.Yi, cfg.Yf, cfg.Ny)
make_omega_2d_input(cfg::EchoConfig) = (t,y) -> pulse_2d_sum(t, y, cfg.pulses)