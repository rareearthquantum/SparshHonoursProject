Base.@kwdef struct EchoConfig
    Nd::Int = 200
    d_width::Float64 = 100.0

    Nt::Int = 500
    Ti::Float64 = 0.0
    Tf::Float64 = 10.0

    Nz::Int = 200
    Zi::Float64 = 0.0
    Zf::Float64 = 10.0

    alpha::Float64 = 1.0e2
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
    twidth = time_window(cfg)
    return t -> pulse(t, twidth/10, twidth/50, pi/2) + pulse(t, 4twidth/10, twidth/100, pi)
end
