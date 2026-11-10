# Input-system roadmap

The project already has validated `EchoConfig` inputs, six Gaussian sequence
presets, Gaussian and top-hat helpers, shared 0D/1D/2D solver support, a transverse
window heuristic, and plotting/JLD2 helpers. Build on these in small steps.
The milestones below are planned work, not implemented interfaces.

1. **Sequence selection.** Add a small configuration builder that resolves named
   recipes using the selected time window and transverse width. Preserve explicit
   custom pulses; test preset resolution and override precedence.
2. **Pulse shape selection.** Represent shape separately from pulse parameters
   and connect the existing Gaussian and top-hat evaluators. Document width and
   area conventions; check analytic values and unchanged Gaussian defaults.
3. **Geometry presets.** Provide convenient 0D, 1D, and 2D choices that resolve to
   the existing axis/count conventions. Test singleton axes, conflicting settings,
   and compatibility with direct `EchoConfig` construction.
4. **Physical inputs and scales.** Define physical inputs separately from solver
   inputs, document the scaling equations, and add explicit nondimensionalisation
   with inverse conversions. Retain scale metadata; test round trips and
   equivalence with manually scaled small calculations.
5. **Automatic time windows.** Infer optional bounds from pulse support/tails and
   expected echo times. Document padding assumptions and preserve explicit bounds;
   test shifted pulses, echo coverage, and truncation reporting.
6. **Automatic grids and spatial windows.** Add optional resolution/domain
   recommendations based on pulse widths, detunings, and diffraction. Reuse the
   existing transverse heuristic where appropriate. Test rules and overrides;
   distinguish recommendations from scientific convergence guarantees.
7. **Resource preview.** Estimate retained arrays, working memory, and output
   sizes before allocation. Report resolved settings without running the solver;
   validate estimates on tiny cases and state runtime-estimate uncertainty.
8. **Configurable plotting and saving.** Add independent output options for
   figures, numerical data, and directories; capture elapsed time and return the
   result. Make research classification optional for general runs. Reuse Plots
   and JLD2, save physical/scaling metadata when available, and test temporary-file
   round trips and disabled outputs.

Each milestone should preserve the existing equations and numerical behaviour,
include focused tests and a short usage example, and honour explicit user settings.
Configuration and saving must use standard Julia utilities and existing
dependencies without DrWatson or another project-management framework. No large
simulation is needed to implement or check the input layer.
