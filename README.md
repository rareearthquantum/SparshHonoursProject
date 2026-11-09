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

The script uses `EchoConfig()` from `src/config_hardcoded.jl`. Its default is a
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
include("src/includes.jl")
cfg = EchoConfig(Nt=128, Nd=32, Nz=8, Ny=0, alpha=100.0, beta=0.0)
result = run_2d_propagation(cfg)
size(result.Omega)  # (time, longitudinal position, transverse position)
size(result.P)      # (time, transverse position)
```

This is an execution example, not a converged scientific calculation.
`run_2d_propagation` returns data without saving or plotting it.
`run_2d_simulation(cfg)` runs the solver and selects the existing research plots.
Include `src/includes.jl` once per Julia session.

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

Presets live in `src/pulse_presets.jl`. Available presets are
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

## Results and source files

`result.Omega[t,z,y]` is the complex field expressed as a Rabi frequency.
`result.P[t,y]` is the detuning-averaged coherence representing polarisation at
the final longitudinal plane by default. With `compute_final_polarisation=false`,
it is from the preceding propagation step; local runs always evaluate it.
Coordinate arrays are `detunings`, `t_vec`, `z_vec`, and `y_vec`.

Within the solver, the descriptive names `rabi_frequency`, `polarisation`, and
`atomic_state` are used. Existing result fields and JLD2 keys remain `Omega` and
`P` so existing scripts and data files keep the same interface.

- `src/config_hardcoded.jl`: configuration defaults and validation.
- `src/input_pulse_methods.jl`: Gaussian and top-hat pulse functions.
- `src/pulse_presets.jl`: example pulse-sequence recipes.
- `src/grids.jl`: coordinate grids and input-field construction.
- `src/integrators.jl`: RK4 atomic integration and AB2 field stepping.
- `src/simulation.jl`: propagation loop and result arrays.
- `src/progress.jl`: progress reporting and elapsed-time formatting.
- `src/plotting/helpers.jl`: shared styles, field transformations, and plot helpers.
- `src/plotting/field_1d.jl`, `field_2d.jl`: field plots for each spatial model.
- `src/plotting/atomic_response.jl`: polarisation and Bloch-component plots.
- `src/plotting/echo_efficiency.jl`: echo-efficiency plots.
- `src/plotting/panels.jl`: composition of plots into multi-panel figures.
- `src/plotting/fitting.jl`: fitted comparisons for SIT, attenuation, and diffraction.
- `src/plotting/selection.jl`: selection and saving of plots for a computed result.
- `src/output/data.jl`: JLD2 result saving.
- `src/output/figures.jl`: static PDF saving.
- `src/output/filenames.jl`: parameter labels shared by output filenames.
- `src/workflow/classification.jl`: experiment classification from its parameters.
- `src/workflow/run_simulation.jl`: running the solver and coordinating figure output.
- `scripts/run_for_plots.jl`: manually selected research examples; inspect before running.
- `notebooks/old/`: historical exploratory notebooks.

`src/includes.jl` remains the main entry point. The existing `src/plotting.jl`,
`src/plot_fitting.jl`, `src/saving.jl`, and `src/sim_run_n_plot.jl` files now load
their corresponding components, preserving existing include paths and function
names. Edit the files in the subfolders when changing their implementations.

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
