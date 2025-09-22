<!--
This file is generated/updated to help AI coding agents (Copilot-style) become productive
in the ShoeOptSetupTime repository. Keep instructions concise and specific to discoverable
patterns in the codebase.
-->

# Copilot instructions for ShoeOptSetupTime

Be brief and precise. Focus on the Julia + DrWatson project layout and optimization solver
patterns used here. Prefer minimal, incremental changes and surface-level tests. When in
doubt, reference the exact files below.

Key files and where to start
- `scripts/run_sims.jl` — entry point used to run experiments. It activates the DrWatson
  project and loads settings from `data/settings/*.jl`.
- `scripts/simulated_annealing.jl` — contains heuristic/search logic referenced by runs.
- `data/settings/H_O2_33.jl` — example settings file that defines problem instances (Pg, Nit,
  temperature schedule, Tl, Gl, etc.). Use it to create reproducible test cases.
- `src/` — small library code and helpers (e.g., `loggers.jl`). Inspect when changing
  logging or experiment wiring.

Big-picture architecture and data flow
- This repository implements an optimization experiment runner (Julia scripts) that:
  1. Activates the project with DrWatson (`@quickactivate "ShoeOptSetupTime"`).
  2. Loads a settings file from `data/settings/` that defines an instance and algorithm
     parameters.
  3. Calls into heuristics/solvers (simulated annealing and a mathematical-programming
     solver) to produce schedules.
  4. Writes results and logs to `data/sims/` (info/debug/warn log files are used by the
     current workflows).

Important patterns and conventions
- Scripts expect local project activation via DrWatson — preserve `using DrWatson` and
  `@quickactivate` calls at the top of scripts unless intentionally changing environment
  activation.
- Settings files are plain Julia scripts that set global variables (e.g., `Pg`, `Nit`,
  `T0`, `Tl`) rather than a config parser. When adding new settings, use the same style
  and document defaults in `README.md`.
- Solver integration: the project currently targets Gurobi (Gurobi.jl). Solver parameters
  such as `Tl` (time limit) are set in settings and propagated to solver call sites.
  Any change to a solver should update both `README.md` and the places that construct
  solver options.

Developer workflows (how to run, test, debug)
- Reproduce locally:
  1. Open Julia REPL
  2. `using Pkg; Pkg.activate("path/to/project"); Pkg.instantiate()`
 3. Run: `julia scripts/run_sims.jl` (uses the settings file indicated inside the script)
- Quick edits and smoke tests: modify a small settings file in `data/settings/` to use
  tiny instances (reduce `Nit`, `Gl`) so runs complete quickly.
- Logs: check `data/sims/` for generated `*_info.log.txt` and other log files when
  reproducing results.

External dependencies and integration points
- Gurobi: the code expects an active Gurobi license by default. If Gurobi is unavailable,
  the project can be switched to HiGHS (`HiGHS.jl`) by adding the package and
  modifying solver construction. Update `README.md` and the solver setup accordingly.
- No networked services or databases are required; experiments are file-based.

Code-change recommendations for AI agents
- Prefer small, testable edits: change one script or helper, run a short example, and
  confirm logs are produced.
- When introducing new CLI flags or settings, add them to an existing settings file
  in `data/settings/` and document usage in `README.md`.
- When modifying solver calls, search for `Gurobi` or `gurobi` strings across the repo
  and change all occurrences to keep behavior consistent.

Examples (concrete suggestions)
- To add HiGHS as an alternative solver, update docs and add a guarded switch in the
  code that constructs the optimizer, for example:

  - Check settings for a `SOLVER` variable (e.g., `SOLVER = :Gurobi` or `:HiGHS`).
  - Use `HiGHS.Optimizer()` when `SOLVER == :HiGHS` and pass solver options similar to
    existing Gurobi options.

- To run a fast smoke test in CI or locally, create a minimal settings file with:
  `Pg = 1; Nit = 1; Gl = 10; Tl = 1` and call `scripts/run_sims.jl`.

What not to change without review
- Do not remove DrWatson activation lines (`@quickactivate`). They are relied upon to
  set the project environment and local paths.
- Avoid changing log file naming or location without updating the README and any log
  readers; tests and users expect outputs in `data/sims/`.

If anything in these instructions is unclear or you need more code examples (e.g., the
exact solver call sites), tell me which area you want expanded and I'll update this
file with concrete code snippets and the exact files to edit.
