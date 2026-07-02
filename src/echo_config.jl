# Configuration and grid helpers.

Base.@kwdef struct EchoConfig
    # Detuning grid
    Nd::Int = 100
    d_width::Float64 = 10.0

    # Time grid
    Nt::Int = 1000
    Ti::Float64 = 0.0
    Tf::Float64 = 10.0

    # Propagation grid
    Nz::Int = 100
    Zi::Float64 = 0.0
    Zf::Float64 = 10.0

    # Field propagation strength
    alpha::Float64 = 1.0
end

function time_window(cfg::EchoConfig)
    return cfg.Tf - cfg.Ti
end

function make_detunings(cfg::EchoConfig)
    return LinRange(-cfg.d_width/2, cfg.d_width/2, cfg.Nd)
end

function make_time_grid(cfg::EchoConfig)
    return LinRange(cfg.Ti, cfg.Tf, cfg.Nt)
end

function make_z_grid(cfg::EchoConfig)
    return LinRange(cfg.Zi, cfg.Zf, cfg.Nz)
end

function make_default_omega_input(cfg::EchoConfig)
    # This expects `pulse` to have been defined already, e.g. by:
    #     include(srcdir("input_pulse_methods.jl"))
    twidth = time_window(cfg)
    return t -> pulse(t, 3twidth/10, twidth/10, 4pi)
end
