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

function default_soliton_2d_pulses(Ti::Real, Tf::Real; y_pulse_width::Real)
    duration = Tf - Ti
    return [
        (PulseParams(Ti + 1duration/10, duration/100, 2pi), PulseParams(0.0, 1.0, 1.0))
    ]
end

function default_weak_2d_pulses(Ti::Real, Tf::Real; y_pulse_width::Real)
    duration = Tf - Ti
    return [
        (PulseParams(Ti + duration/2, duration/10, pi/8), PulseParams(0.0, 1.0, 1.0))
    ]
end

function default_pi_2d_pulses(Ti::Real, Tf::Real; y_pulse_width::Real)
    duration = Tf - Ti
    return [
        (PulseParams(Ti + 3duration/10, duration/50, pi), PulseParams(0.0, 1.0, 1.0))
    ]
end

Base.@kwdef struct EchoConfig
    Nt::Int = 256 * 2
    Ti::Float64 = 0.0
    Tf::Float64 = 1.0

    d_width::Float64 = 1000
    Nd::Int = 64 * 2 * 2

    Nz::Int = 64 * 2 
    Zi::Float64 = 0.0
    Zf::Float64 = 1.0

    alpha::Float64 = 3000.0
    beta::Float64 = 0.0

    Ny::Int = 64 * 2
    y_pulse_width::Float64 = 1.0
    Yi::Float64 = -4*y_pulse_width
    Yf::Float64 = 4*y_pulse_width

    pulses::Vector{NTuple{2,PulseParams}} = default_echo_2d_pulses(Ti, Tf; y_pulse_width)

    function EchoConfig(
        Nt, Ti, Tf,
        d_width, Nd,
        Nz, Zi, Zf,
        alpha, beta,
        Ny, y_pulse_width, Yi, Yf,
        pulses,
    )
        Nt = Int(Nt)
        Nd = Int(Nd)
        Nz = Int(Nz)
        Ny = Int(Ny)

        Ti = Float64(Ti)
        Tf = Float64(Tf)
        d_width = Float64(d_width)
        Zi = Float64(Zi)
        Zf = Float64(Zf)
        alpha = Float64(alpha)
        beta = Float64(beta)
        y_pulse_width = Float64(y_pulse_width)
        Yi = Float64(Yi)
        Yf = Float64(Yf)

        Nt >= 1 || throw(ArgumentError("Nt must be positive"))
        Nd >= 1 || throw(ArgumentError("Nd must be positive"))
        Nz >= 1 || throw(ArgumentError("Nz must be positive"))
        Ny >= 1 || throw(ArgumentError("Ny must be positive"))

        Ti < Tf || throw(ArgumentError("Ti must be less than Tf"))
        Zi < Zf || throw(ArgumentError("Zi must be less than Zf"))
        Ny == 1 || Yi < Yf ||
            throw(ArgumentError("Yi must be less than Yf when Ny > 1"))

        d_width > 0 ||
            throw(ArgumentError("d_width must be positive"))
        y_pulse_width > 0 ||
            throw(ArgumentError("y_pulse_width must be positive"))

        Ny == 1 || Nz > 1 ||
            throw(ArgumentError("Nz must be greater than 1 when Ny > 1"))

        return new(
            Nt, Ti, Tf,
            d_width, Nd,
            Nz, Zi, Zf,
            alpha, beta,
            Ny, y_pulse_width, Yi, Yf,
            pulses,
        )
    end
end

make_detunings(cfg::EchoConfig) = LinRange(-cfg.d_width/2, cfg.d_width/2, cfg.Nd)
make_time_grid(cfg::EchoConfig) = LinRange(cfg.Ti, cfg.Tf, cfg.Nt)
make_z_grid(cfg::EchoConfig) = (cfg.Nz==1) ? [0.0] : LinRange(cfg.Zi, cfg.Zf, cfg.Nz)
make_y_grid(cfg::EchoConfig) = (cfg.Ny==1) ? [0.0] : LinRange(cfg.Yi, cfg.Yf, cfg.Ny)
make_ky_grid(cfg::EchoConfig) = (cfg.Ny==1) ? [0.0] : 2pi .* fftfreq(cfg.Ny, (cfg.Ny - 1) / (cfg.Yf - cfg.Yi))
make_omega_2d_input(cfg::EchoConfig) = (t, y) -> pulse_2d_sum(t, y, cfg.pulses)