# SparshHonoursProject

## Introduction

This repository contains code for Sparsh Chandra's PHSI490 research project,
supervised by Associate Professor Jevon Longdell at the University of Otago,
for the partial completion of a Bachelor of Science with Honours in Physics.
The simulations model a detuned two-level atomic ensemble driven by optical
pulses, with longitudinal propagation and optional transverse diffraction.
Examples include photon echoes and self-induced transparency (SIT).

## Setup

Install Git and Julia. The checked-in package environment uses **Julia 1.12.6**.
From a terminal:

```sh
git clone https://github.com/rareearthquantum/SparshHonoursProject.git
cd SparshHonoursProject
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The initial package download and compilation can take several minutes. Keep
`Project.toml` and `Manifest.toml` together when reproducing a run. Run the
commands below from the repository root.

## Running

```sh
julia --project=. scripts/run_simulation.jl
```

The script uses `EchoConfig()` from `src/config/echo_config.jl`. Its default is a
weak Gaussian pulse with diffraction and no atom–field coupling (`alpha=0`,
`beta=1`). Figures are written under `plots/results/testing/finalstretch/`.
**The script currently saves plots only**; numerical-data saving is disabled in
`src/workflow/run_simulation.jl`.

For a Linux server without a display:

```sh
GKSwstype=100 julia --project=. scripts/run_simulation.jl
```

To try a small calculation interactively, start `julia --project=.` and enter:

```julia
using Dates
include("src/main_include.jl")
cfg = EchoConfig(Nt=128, Nd=32, Nz=8, Ny=0, alpha=100.0, beta=0.0)
result = run_2d_propagation(cfg)
size(result.Omega)  # (time, longitudinal position, transverse position)
size(result.P)      # (time, transverse position)
```

This is an execution example, not a converged scientific calculation.
`run_2d_propagation` returns data without saving or plotting it.
`run_2d_simulation(cfg)` runs the solver and selects the existing research plots.
Include `src/main_include.jl` once per Julia session.

Progress reports actual completed work at roughly 10% intervals, including the
final polarisation calculation. Use `run_2d_propagation(cfg; show_progress=false)`
to suppress progress messages.

## Configuration

Set parameters using `EchoConfig` keyword arguments in your run script:

| Parameter | Meaning |
| --- | --- |
| `Nt`, `Ti`, `Tf` | Number of time samples and time interval; `Nt >= 2`, `Ti < Tf` |
| `Nd`, `d_width` | Number of equally weighted detunings and their full interval width, centred on zero |
| `Nz`, `Zi`, `Zf` | Longitudinal sample count and interval |
| `Ny`, `Yi`, `Yf` | Transverse sample count and interval |
| `alpha` | Atom–field coupling coefficient |
| `beta` | Transverse diffraction coefficient |
| `y_pulse_width` | Transverse Gaussian standard deviation used by presets and default transverse bounds |
| `pulses` | Vector of `(temporal, transverse)` `PulseParams` pairs |

`Nd`, `Nz`, and `Ny` may each be **0 or 1 to select a single point**. Both values
are stored as 1 in the configuration. A single detuning is zero; a single spatial
coordinate is zero. Negative counts are rejected. Time evolution needs at least
two samples, so `Nt=0` and `Nt=1` are rejected.

| Spatial model | Settings |
| --- | --- |
| Local atomic response (0D) | `Nz=0` or `1`, and `Ny=0` or `1` |
| Longitudinal propagation (1D) | `Nz>1`, and `Ny=0` or `1` |
| Longitudinal and transverse grid (2D) | `Nz>1`, `Ny>1` |

A transverse grid with more than one point requires longitudinal propagation.
Spatial interval endpoints must increase when that axis has multiple points;
single-point spatial axes ignore their interval. Detuning width must be positive
for multiple detunings, and may be zero for a single detuning. Pulse widths and
`y_pulse_width` must be positive. Scalar and pulse parameters must be finite.

Presets live in `src/pulses/presets.jl`. Available presets are
`default_weak_2d_pulses`, `default_pi_2d_pulses`,
`default_soliton_2d_pulses`, `default_4pi_2d_pulses`, `default_6pi_2d_pulses`, and
`default_echo_2d_pulses`. All accept `Ti`, `Tf`, and `y_pulse_width` keywords.
When supplying a preset explicitly, pass the same time window and transverse
width to it and the configuration:

```julia
cfg = EchoConfig(
    Tf=2.0, y_pulse_width=0.5, alpha=3000.0, beta=0.0, Ny=1,
    pulses=default_echo_2d_pulses(Tf=2.0, y_pulse_width=0.5),
)
```

Temporal Gaussian `width` is a standard deviation and `area` is the integral over
the full real line; a finite time window may truncate it. Transverse Gaussians
have unit peak, use `width` as their standard deviation, and ignore their `area`
field. Explicit pulses keep their supplied parameters when other config values
change.

The current presets are example pulse recipes, not full experiment definitions.
For a future interface, a single Gaussian-sequence builder with explicit area,
centre, and width would cover the single-pulse cases; a separate echo builder
would describe the two-pulse sequence. A configuration keyword such as
`pulse_preset=:echo` could then resolve the recipe using that configuration's
`Ti`, `Tf`, and `y_pulse_width`, avoiding the duplicated arguments above.
That keyword is a proposed interface, not an implemented option. Explicit
`pulses` should remain available for custom sequences.

The existing presets intentionally use different timings and widths: changing
only the area of a weak pulse does not reproduce the SIT or 4pi examples. Names
such as `soliton` also describe an intended experiment, rather than guaranteeing
that a pulse will propagate as a soliton for every grid and coupling strength.

Use a consistent numerical unit system. Time and detuning have reciprocal units,
with detuning interpreted as angular frequency because it appears in
`cis(detuning*t)`. `Omega` is a Rabi frequency with inverse-time units. From the
implemented equations, `alpha` has inverse-time per longitudinal-length units,
and `beta` has transverse-length squared per longitudinal-length units.
All parameters are plain numeric values; no unit package or automatic conversion
is used by the simulation. Record the chosen scales with shared results.

## Constructing user-facing inputs

`SimulationConfig`, `PulseSpec`, `MediumSpec`, and `GridSpec` are standalone
input records. Load them explicitly; the current entry points continue to use
`EchoConfig`. They support construction, field access, and standard Julia display
only. There is no conversion to solver inputs, unit conversion, grid generation,
preset expansion, or physical/cross-field validation yet.

For example, from a Julia session in the project root:

```julia
include("src/config/simulation_config.jl")

# Here we choose seconds, metres, and angular frequencies in rad/s.
cfg = SimulationConfig(
    geometry=:one_d,
    pulses=[PulseSpec(center=2e-6, width=0.2e-6, area=pi/2, y_width=0.5e-3)],
    medium=MediumSpec(d_width=2pi*1e6, alpha=2e8, beta=0.0),
    grid=GridSpec(Nt=64, Ti=0.0, Tf=5e-6, Nd=16,
        Nz=8, Zi=0.0, Zf=0.01, Ny=1, Yi=-0.002, Yf=0.002),
)

@show cfg
@show cfg.pulses[1].width cfg.medium.d_width cfg.grid.Nt
```

The example coefficients illustrate input syntax rather than a calibrated
material. Numerical defaults use the existing project's familiar values, without
assuming SI units. Gaussian temporal widths are standard deviations, temporal
areas are full-pulse integrals of Rabi frequency, and transverse profiles are
described as unit-peak Gaussians. Geometry and shape symbols are descriptive
labels; they do not apply presets or rewrite settings. Default pulses and bounds
are fixed values rather than functions of other fields. Each record's docstring
describes its fields and units.

Construction-only checks can run without loading the solvers:

```sh
julia --project=. test/configuration_construction.jl
```

## Results and source files

`result.Omega[t,z,y]` is the complex field expressed as a Rabi frequency.
`result.P[t,y]` is the detuning-averaged coherence representing polarisation at
the final longitudinal plane by default. With `compute_final_polarisation=false`,
it is from the preceding propagation step; local runs always evaluate it.
Coordinate arrays are `detunings`, `t_vec`, `z_vec`, and `y_vec`.

Within the solver, the descriptive names `rabi_frequency`, `polarisation`, and
`atomic_state` are used. Existing result fields and JLD2 keys remain `Omega` and
`P` so existing scripts and data files keep the same interface.

`src/main_include.jl` is the main entry point. Files ending in `_include.jl`
only load the files in their component, in dependency order.

```text
src/
├── main_include.jl
├── config/
│   ├── echo_config.jl       # Existing solver configuration and validation
│   └── simulation_config.jl # Standalone user-facing input records
├── pulses/
│   ├── parameters.jl        # PulseParams and its constructor
│   ├── gaussian.jl          # Gaussian evaluations
│   ├── top_hat.jl           # Top-hat evaluations
│   ├── sequences.jl         # Pulse sums and separable 2D inputs
│   └── presets.jl           # Example experiment recipes
├── integrators/
│   ├── adams_bashforth.jl   # AB2 cache and stepping
│   └── atomic_rk4.jl        # RK4 cache and atomic-state integration
├── simulation/
│   ├── grids.jl             # Coordinate grids and input-field construction
│   ├── atomic_response.jl   # Coherence phases and ensemble polarisation
│   ├── field.jl             # Spectral field derivative
│   ├── propagation.jl       # Propagation loop and result arrays
│   └── progress.jl          # Completion reporting and timing
├── plotting/
│   ├── style.jl            # Shared sizes, margins, colours, and labels
│   ├── field_data.jl       # Field transformations and threshold bounds
│   ├── layout.jl           # Panel labels and figure composition
│   ├── primitives.jl       # Shared axes, heatmaps, and temporal curves
│   ├── field_1d.jl         # Longitudinal field plots
│   ├── field_2d.jl         # Transverse field plots
│   ├── atomic_response.jl  # Polarisation and Bloch-component plots
│   ├── echo_efficiency.jl  # Echo-efficiency estimates and plots
│   ├── panels.jl           # Experiment-specific figure combinations
│   └── fitting/            # Fit models and SIT/attenuation/diffraction plots
├── output/                 # Filename labels, JLD2 data, and PDF output
└── workflow/               # Experiment classification, plot selection, and execution
```

Each multi-file component has a corresponding `_include.jl` loader, omitted
from the tree above for readability. The tests load the solver and output
components without loading the plotting dependencies. Run scripts still use
`src/main_include.jl`; experiment configuration belongs in those scripts.
`scripts/run_for_plots.jl` contains manually selected research examples, so
inspect it before running.

The root `Project.toml` declares Julia dependencies, and `Manifest.toml` locks
versions for reproducible runs. `.gitignore` excludes generated output and
temporary files; `.gitattributes` controls text handling and line endings.
`.harper-dictionary.txt` supplies project-specific words to the Harper spellchecker.
Git keeps local repository metadata in the hidden `.git/` directory.

Generated files in `data/` and `plots/` are ignored by Git. For shared results,
retain the run parameters, unit conventions, Julia environment, code revision,
and convergence checks. The transverse solver loop is currently serial.

## Cluster jobs and checks

`scripts/run_echo_2d.slurm` invokes `scripts/run_simulation.jl` with that script's
configuration, despite the historical Slurm filename. Adapt the partition,
resources, and notification settings to your cluster. Julia must be available on
the compute node and the environment must already be instantiated.

Create the log directory **before submission**, because Slurm opens the log
before running the job script:

```sh
mkdir -p data/slurm
sbatch scripts/run_echo_2d.slurm
```

Run the small configuration and solver regression checks with:

```sh
julia --project=. test/runtests.jl
```

These checks do not establish scientific convergence. Check time, longitudinal,
detuning, and transverse resolution and domain size for each experiment.
