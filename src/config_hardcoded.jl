y_width(y_pulse_width, beta, Z_length) = 6*y_pulse_width*sqrt(1+(2*beta*Z_length/y_pulse_width^2)^2)

"""
    EchoConfig(; kwargs...)

Simulation parameters. `Nd`, `Nz`, and `Ny` values of 0 or 1 select a single
point and are stored as 1. Time evolution requires `Nt >= 2`.
Single detuning and spatial coordinates are zero. See README.md for units,
pulse conventions, and the supported spatial grids.
"""
Base.@kwdef struct EchoConfig
    Nt::Int = 256
    Ti::Float64 = 0.0
    Tf::Float64 = 1.0

    d_width::Float64 = 1000
    Nd::Int = 128

    Nz::Int = 64
    Zi::Float64 = 0.0
    Zf::Float64 = 1.0

    alpha::Float64 = 0.0
    beta::Float64 = 1.0

    Ny::Int = 128
    y_pulse_width::Float64 = 1.0
    Yi::Float64 = -y_width(y_pulse_width, beta, Zf-Zi)/2
    Yf::Float64 = y_width(y_pulse_width, beta, Zf-Zi)/2

    pulses::Vector{NTuple{2,PulseParams}} = default_weak_2d_pulses(; Ti, Tf, y_pulse_width)


    """
    This following constructor just checks the input config and makes sure
    the values are valid
    """

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

        Nt >= 2 || throw(ArgumentError("Nt must be at least 2 to define a time step"))
        Nd >= 0 || throw(ArgumentError("Nd must be non-negative"))
        Nz >= 0 || throw(ArgumentError("Nz must be non-negative"))
        Ny >= 0 || throw(ArgumentError("Ny must be non-negative"))
        Nd = max(1, Nd)
        Nz = max(1, Nz)
        Ny = max(1, Ny)

        all(isfinite, (Ti, Tf, d_width, Zi, Zf, alpha, beta, y_pulse_width, Yi, Yf)) ||
            throw(ArgumentError("configuration values must be finite"))
        for pair in pulses, params in pair
            all(isfinite, (params.center, params.width, params.area)) ||
                throw(ArgumentError("pulse parameters must be finite"))
            params.width > 0 || throw(ArgumentError("pulse widths must be positive"))
        end

        Ti < Tf || throw(ArgumentError("Ti must be less than Tf"))
        Nz == 1 || Zi < Zf ||
            throw(ArgumentError("Zi must be less than Zf when Nz > 1"))
        Ny == 1 || Yi < Yf ||
            throw(ArgumentError("Yi must be less than Yf when Ny > 1"))

        (Nd == 1 ? d_width >= 0 : d_width > 0) ||
            throw(ArgumentError("d_width must be positive for Nd > 1, or non-negative for a single detuning"))
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
