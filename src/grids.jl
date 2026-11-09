using FFTW: fftfreq

# Grid helpers and conversion from a configuration to numerical coordinates.
const ND_TO_NT_RATIO = 0.2

detuning_width(Nt, Ti, Tf) = (Nt - 1) / (Tf - Ti) # 1/dt
detuning_count(Nt) = ceil(Int, ND_TO_NT_RATIO * Nt)

make_detunings(cfg::EchoConfig) = cfg.Nd == 1 ? [0.0] : LinRange(-cfg.d_width/2, cfg.d_width/2, cfg.Nd)
make_time_grid(cfg::EchoConfig) = LinRange(cfg.Ti, cfg.Tf, cfg.Nt)
make_z_grid(cfg::EchoConfig) = (cfg.Nz==1) ? [0.0] : LinRange(cfg.Zi, cfg.Zf, cfg.Nz)
make_y_grid(cfg::EchoConfig) = (cfg.Ny==1) ? [0.0] : LinRange(cfg.Yi, cfg.Yf, cfg.Ny)
make_ky_grid(cfg::EchoConfig) = (cfg.Ny==1) ? [0.0] : 2pi .* fftfreq(cfg.Ny, (cfg.Ny - 1) / (cfg.Yf - cfg.Yi))
make_omega_2d_input(cfg::EchoConfig) = (t, y) -> pulse_2d_sum(t, y, cfg.pulses)
