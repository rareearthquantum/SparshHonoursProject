using FFTW: fftfreq

# Grid helpers and conversion from a configuration to numerical coordinates.
const ND_TO_NT_RATIO = 0.2

"""
    detuning_width(Nt, Ti, Tf)

Return the inverse time spacing `(Nt - 1) / (Tf - Ti)` for an inclusive
uniform time grid. Requires `Nt >= 2` and `Tf > Ti`.
"""
detuning_width(Nt, Ti, Tf) = (Nt - 1) / (Tf - Ti) # 1/dt
"""
    detuning_count(Nt)

Return the default detuning sample count, rounding `ND_TO_NT_RATIO * Nt`
up to an integer.
"""
detuning_count(Nt) = ceil(Int, ND_TO_NT_RATIO * Nt)

"""
    make_detunings(cfg::EchoConfig)

Return `cfg.Nd` equally weighted detunings spanning
`[-cfg.d_width/2, cfg.d_width/2]`, or `[0.0]` for a single detuning.
"""
make_detunings(cfg::EchoConfig) = cfg.Nd == 1 ? [0.0] : LinRange(-cfg.d_width/2, cfg.d_width/2, cfg.Nd)
"""
    make_time_grid(cfg::EchoConfig)

Return an inclusive uniform time range from `cfg.Ti` to `cfg.Tf`
with `cfg.Nt` samples.
"""
make_time_grid(cfg::EchoConfig) = LinRange(cfg.Ti, cfg.Tf, cfg.Nt)
"""
    make_z_grid(cfg::EchoConfig)

Return an inclusive uniform longitudinal range with `cfg.Nz` samples,
or `[0.0]` when the longitudinal axis has one point.
"""
make_z_grid(cfg::EchoConfig) = (cfg.Nz==1) ? [0.0] : LinRange(cfg.Zi, cfg.Zf, cfg.Nz)
"""
    make_y_grid(cfg::EchoConfig)

Return an inclusive uniform transverse range with `cfg.Ny` samples,
or `[0.0]` when the transverse axis has one point.
"""
make_y_grid(cfg::EchoConfig) = (cfg.Ny==1) ? [0.0] : LinRange(cfg.Yi, cfg.Yf, cfg.Ny)
"""
    make_ky_grid(cfg::EchoConfig)

Return transverse angular wavenumbers in FFT order, using the spacing
of the configured transverse grid. Return `[0.0]` for a single point.
"""
make_ky_grid(cfg::EchoConfig) = (cfg.Ny==1) ? [0.0] : 2pi .* fftfreq(cfg.Ny, (cfg.Ny - 1) / (cfg.Yf - cfg.Yi))
"""
    make_omega_2d_input(cfg::EchoConfig)

Return a callable `(t, y) -> field` that evaluates the Gaussian pulse
sequence stored in `cfg.pulses`.
"""
make_omega_2d_input(cfg::EchoConfig) = (t, y) -> pulse_2d_sum(t, y, cfg.pulses)
