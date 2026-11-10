# Standalone user-input records. Load explicitly; existing entry points continue
# to use EchoConfig. No solver, pulse evaluators, or other dependencies are loaded.

"""
    PulseSpec(; shape=:gaussian, center=0.5, width=0.1, area=pi/8,
                y_center=0.0, y_width=1.0)

Describe one temporal pulse and its unit-peak Gaussian transverse profile.
`shape` labels the temporal shape (for example `:gaussian` or `:top_hat`).
`center` and `width` use the caller's time units; `area` is the dimensionless
Rabi-frequency integral over the full pulse. Gaussian `width` means standard
deviation; top-hat `width` means full support width. `y_center` and `y_width`
use transverse length units, with `y_width` the Gaussian standard deviation.

Defaults describe the existing weak-pulse recipe on a 0-to-1 time window,
but are fixed values: changing a window does not move or resize this pulse.
This record stores parameters only; it does not evaluate or expand a pulse.
"""
Base.@kwdef struct PulseSpec
    shape::Symbol = :gaussian
    center::Float64 = 0.5
    width::Float64 = 0.1
    area::Float64 = pi/8
    y_center::Float64 = 0.0
    y_width::Float64 = 1.0
end

"""
    MediumSpec(; d_width=1000.0, alpha=0.0, beta=1.0)

Describe a medium using the existing effective-coefficient naming conventions.
`d_width` is the full angular-detuning interval width, with inverse-time units.
`alpha` has units of inverse time per longitudinal length; `beta` has units of
transverse length squared per longitudinal length. Defaults follow EchoConfig.

These are user-supplied dimensionful coefficients in a consistent unit system,
not automatically nondimensionalised solver values. No material constants are
inferred and no conversions or detuning samples are computed.
"""
Base.@kwdef struct MediumSpec
    d_width::Float64 = 1000.0
    alpha::Float64 = 0.0
    beta::Float64 = 1.0
end

"""
    GridSpec(; Nt=256, Ti=0.0, Tf=1.0, Nd=128,
               Nz=64, Zi=0.0, Zf=1.0, Ny=128, Yi=-3.0, Yf=3.0)

Store requested sample counts and explicit time, longitudinal, and transverse
bounds, using the existing EchoConfig names. Bounds use the caller's respective
time and length units. `Nd` is the detuning sample count; its interval width is
stored in MediumSpec. The default transverse bounds are fixed at -3 and 3;
they are not calculated from pulse widths or diffraction.

No coordinate arrays, automatic bounds, or resolutions are generated. Counts
are stored as supplied, including zero; there is no singleton-axis normalisation.
"""
Base.@kwdef struct GridSpec
    Nt::Int = 256
    Ti::Float64 = 0.0
    Tf::Float64 = 1.0
    Nd::Int = 128
    Nz::Int = 64
    Zi::Float64 = 0.0
    Zf::Float64 = 1.0
    Ny::Int = 128
    Yi::Float64 = -3.0
    Yf::Float64 = 3.0
end

"""
    SimulationConfig(; geometry=:two_d, pulses=[PulseSpec()],
                       medium=MediumSpec(), grid=GridSpec())

Group user-facing inputs for construction and inspection. `geometry` is a label
such as `:zero_d`, `:one_d`, or `:two_d`; it does not change counts, coefficients,
or pulses. `pulses` is an explicit sequence, with one fixed weak Gaussian by
default; it may also be empty. Each default configuration gets its own vector.

All numerical values use a caller-chosen consistent unit system; SI is not
assumed. These data-only records do not validate physical or cross-field
consistency, convert units, resolve presets, or convert to EchoConfig. Existing
solvers and entry points continue to accept EchoConfig only. Standard Julia
field access and display can be used to inspect the records.
"""
Base.@kwdef struct SimulationConfig
    geometry::Symbol = :two_d
    pulses::Vector{PulseSpec} = [PulseSpec()]
    medium::MediumSpec = MediumSpec()
    grid::GridSpec = GridSpec()
end
