# SparshHonoursProject

## Introduction

This repository contains the code for the PHSI490 research project of Sparsh Chandra, supervised by [Associate Professor Jevon Longdell](https://www.otago.ac.nz/physics/staff/jevonlongdell) at the University of Otago, for the partial completion of a Bachelor of Science with Honours (BSc(Hons)) in Physics.

## Prerequisites

You need to have [Git](https://git-scm.com/) and [Julia](https://julialang.org/) installed on your computer.
This project is fully runnable on either Windows or Linux - macOS has yet to be tested.

## First Install

To (locally) reproduce this project, do the following:
(NOTE: the `instantiate` command might take some time to precompile the packages depending on your internet connection and computer specs; up to 20+ minutes):
```
terminal> cd path/to/where/you/want/to/put/it/
terminal> git clone https://github.com/rareearthquantum/SparshHonoursProject.git
terminal> cd SparshHonoursProject
julia> ]
pkg> activate .
pkg> instantiate
```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

## Running

Now once its set up, to run any script from a fresh terminal session:

```
terminal> cd path/to/SparshHonoursProject
terminal> julia --project=. ./scripts/run_simulation.jl

```

## Configuring Run

In src/config_hardcoded.jl will be the different parameters.

For 0SD (no spatial propagation at all) set Nz and Ny to 1.
For 1SD (no transverse propagation) set Nz>1 and Ny=1.
For 2SD both Nz and Ny > 1.

To change plots, check src/plotting.jl, make sure the right plotting functions are defined or define your own, then call them in src/sim_run_n_plot.jl.