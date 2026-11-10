# Simulation regression checks

Run the full suite from the project root:

```sh
julia --project=. test/runtests.jl
```

`simulation_regressions.jl` is included by that runner. It adds these small cases:

| Case | Grid counts | Reference or invariant |
|---|---|---|
| Local atoms | Nt=129, Nd=1 or 5, Nz=Ny=1 | Exact constant-drive Rabi coherence, including detuning and a shifted time origin |
| Undriven atoms | Nt=17, Nd=5, Nz=Ny=1 | Zero field and coherence remain zero |
| Coupled 1D | Nt=65, Nd=Ny=1, Nz=9 then 17 | Weak-field analytical solution and second-order longitudinal refinement |
| Pure diffraction | Nt=17, Nd=3, Nz=5, Ny=8 | Exact periodic Fourier-mode phases and transverse squared-norm conservation |
| Coupled 2D | Nt=65, Nd=3, Nz=5, Ny=8 | Uniform transverse input agrees with a corresponding 1D calculation |

Tolerances are absolute maximum errors per entry, justified beside the assertions.
Rabi checks allow RK4 truncation error; weak-field checks allow the small nonlinear
correction and spatial discretisation error. Spectral diffraction and dimensional
reduction checks allow floating-point roundoff. References are analytical, with no
saved solver snapshots or new dependencies. The diffraction domain contains eight
distinct periodic samples, avoiding a duplicated endpoint.

All cases belong in the default suite. Julia loading and compilation can dominate
elapsed time. Large echo/SIT runs and broad convergence studies should remain
optional, separately approved work; these small checks do not establish convergence
for research configurations. No plotting, cluster jobs, or result files are needed
by these added cases; the existing saving tests use temporary files.
