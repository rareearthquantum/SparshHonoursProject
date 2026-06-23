# SparshHonoursProject

## Introduction

This repository contains the code and dissertation for the PHSI490 research project of Sparsh Chandra. Supervised by [Associate Professor Jevon Longdell](https://www.otago.ac.nz/physics/staff/jevonlongdell) at the University of Otago, for the partial completion of a Bachelor of Science with Honours (BSc(Hons)) in Physics.

This project aims to model photon echoes in 2 spatial dimensions (2SD). The area theorem tells us that the time integral of a pulse dictates the evolution of the Bloch vectors of each atom due to the pulse. However, this integral is only in time, so it can vary in space. Most laser beams - which will be used to produce photon echoes - are radially Gaussian, meaning the area of the pulse will be different in the transverse direction, affecting the way photon echoes are produced. Hence, having a numeric model and/or simulation would help in the exploration of such systems.

This project is working from the ground up. First starting off with the Maxwell-Bloch equations in input-output theory formalism, and then treating the operators as scalars (to begin with a more classical case) and letting $\sigma_z \approx -1$ (linearising the equations), all still in 1SD. Then as the project progresses, more of the original equations will be reintroduced, while also extending it into 2SD.

2SD, scalars, linear: this is the current case, where we will be working to stabilise the current code and then extend from scalars -> operators. An extension in 3SD may even be possible in a centrosymmetric case.

## Prerequisites

You need to have [Git](https://git-scm.com/) and [Julia](https://julialang.org/) installed on your computer.
This project is fully runnable on either Windows or Linux - macOS has yet to be tested.

### Technical Lingo/Jargon/Words for the Uninitiated

 - Git: A version control system and software, designed to track changes in code or files
 - GitHub: A cloud-based service that allows you to manage, share and collaborate on Git repositories online
 - Repository (repo): The database that stores your project's files and complete history of changes
 - Clone: Copies an existing remote repository from a platform like GitHub to your local computer
 - Julia: A programming language
 - Pkg: Julia's built in package manager which handles installing, updating, and removing packages
 - Packages: Reusable bundle of code which can extend the functionality of Julia by providing new functions and more
 - [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/): A Julia package designed to make the lives of scientists easier by making scientific projects easy to reproduce

## First Install

To (locally) reproduce this project, do the following:
(NOTE: the `instantiate` command might take some time to precompile the packages depending on your internet connection and computer specs; up to 20+ minutes):
```
terminal> cd path/to/where/you/want/to/put/it/
terminal> git clone https://github.com/rareearthquantum/SparshHonoursProject.git
terminal> cd SparshHonoursProject
julia> ]
pkg> add DrWatson
pkg> activate .
pkg> instantiate
```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

## Running

Now once its set up, to run any script from a fresh terminal session:

### Option 1

Nice and simple

```
terminal> cd path/to/SparshHonoursProject
terminal> julia scripts/run.jl
```

And to rerun:

```
terminal> julia scripts/run.jl
```

### Option 2

Access to Julia REPL and Pkg, also successive reruns take less time:

```
terminal> cd path/to/SparshHonoursProject
terminal> julia
julia> include("scripts/run.jl")
```

And to rerun:

```
julia> include("scripts/run.jl")
```

### Prototype configurations

The FFT prototype reads a TOML configuration. With no argument it uses the
dimensionless single-pulse configuration:

```text
julia --project=. scripts/prototype_fft_run.jl
```

Pass the physical SI configuration explicitly with:

```text
julia --project=. scripts/prototype_fft_run.jl scripts/input_params_2_10micron_y_pulses_40micron_seperation.toml
```

TOML files contain data only: expressions such as `10micro`, `128 * 2`, and
Julia tuples are not valid TOML. The physical configuration therefore uses
plain SI values with unit suffixes such as `_m`, `_s`, and `_hz_per_m`.

The numerical kernel is always dimensionless. For reference scales `Z0`, `T0`
and `Y0`, the loader uses

```text
alpha = coupling_hz_per_m * Z0 * T0
beta  = (1 / (2k)) * Z0 / Y0^2
```

and converts all domains and pulse parameters by the same scales. This keeps
units out of the FFT and integration kernels while restoring SI coordinates on
physical-mode plot axes.
