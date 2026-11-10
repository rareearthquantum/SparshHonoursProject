# Coding agent guidance

- Inspect the relevant source and tests before editing; check whether a requested
  feature already exists. See `README.md` and `ROADMAP.md` for context.
- Preserve the existing Maxwell–Bloch equations, conventions, and numerical
  behaviour unless the user explicitly requests a change.
- Prefer small, independently testable changes. Avoid unnecessary dependencies
  and large refactors; reuse existing Julia code and conventions where practical.
- Keep user-facing physical parameters separate from dimensionless solver
  parameters. Make conversions explicit and retain the scales needed to interpret
  results. Do not assume existing plain numerical inputs are in SI units.
- This project does not use DrWatson. Do not introduce it or use its path,
  naming, saving, or experiment-management functions, even if a task mentions it.
  Use standard Julia utilities such as `joinpath`, `@__DIR__`, and `mkpath`, and
  existing dependencies such as JLD2. Keep configuration and saving independent
  of project-management frameworks.
- Never silently change explicit user settings. Presets and automatic choices
  must respect overrides; report invalid or conflicting inputs clearly.
- Explain mathematical assumptions and non-obvious numerical choices, including
  units, pulse width/area conventions, boundary conditions, and grid heuristics.
- Add or update tests for behaviour-changing code and run relevant small checks.
  The existing suite is `julia --project=. test/runtests.jl` from the project root.
  Report what was checked and any checks that were not run.
- Do not run large simulations, parameter sweeps, or cluster jobs without
  explicit user approval. Run scripts may launch substantial calculations.
- Never automatically commit or push changes; do so only when explicitly asked.
