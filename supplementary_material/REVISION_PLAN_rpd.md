# Revision Plan — supplementary_rpd.md (reviewer Q2 & Q3)

Goal: make `supplementary_rpd.md` fully answer reviewer requests #2 (average RPD
per instance for all stochastic methods, BKS = best across all exact+heuristic
methods) and #3 (rigorous statistical analysis: ≥30 runs, normality,
parametric vs nonparametric, post-hoc).

Source of this plan: code review on 2026-06-01.

---

## Critical (must fix — will trigger another review round otherwise)

- [ ] **State the number of independent runs.** Add to the methods note the exact
  number of runs per instance per stochastic method (must be ≥30). This is the
  core of Q3 and is currently missing entirely.

- [ ] **Resolve SSM-SA inconsistency.** SSM-SA is described as stochastic but
  reported as a single value (no std), and `—` is explained as "infeasible for
  repeated runs" (line 18). Decide and document one of:
    - (a) report SSM-SA avg + std over the same run count as the others, OR
    - (b) keep single-run, but state it explicitly, defend why, and reconsider
      whether it belongs in the Friedman/Wilcoxon tests (mixing single-sample
      with 30-run averages is questionable).

- [ ] **Add the constructive (greedy) baseline.** Line 35 mentions a
  "deterministic baseline" but there is no column for it. Add its RPD as a
  reference column/row, since the reviewer explicitly contrasts the greedy
  constructive algorithm with the stochastic ones.

## Internal consistency

- [ ] **Matched instance set for the summary.** Summary table (lines 25–31) mixes
  n=120 (SA/GRASP/GA), n=109 (SSM-SA), n=94 (MILP), so the headline means are not
  comparable. Either add a second summary over the common 109-instance set, or add
  an explicit caveat that means span different subsets.

- [ ] **Normality on the matched set.** Shapiro–Wilk uses n=120 (lines 175–178)
  while Friedman/Wilcoxon use n=109. Re-run normality on the same 109 matched set.

- [ ] **State MILP's role.** Add one line: statistical analysis compares the four
  stochastic/heuristic methods; MILP is the exact reference and is excluded from
  the metaheuristic comparison.

## Smaller fixes

- [ ] Typo line 207: "is available is" → "is available below".
- [ ] Confirm `rpd_boxplot.png` is committed next to the file.
- [ ] Add the reviewer-requested citation: Molina et al., *Swarm and Evolutionary
  Computation* 64, 100888 (DOI: 10.1016/j.swevo.2021.100888), alongside Derrac et al. (2011).
- [ ] State the effect-size thresholds used (Cohen r: 0.1/0.3/0.5) so labels are
  reproducible; re-check SA vs SSM-SA r=0.579 "large" label against them.
- [ ] For rows with std=0 and nonzero mean (e.g. lines 42, 102): confirm these are
  real convergence over the full run count, not single runs (ties to run-count fix).
- [ ] State seeds + run count in methods and confirm source code/seeds are in the
  GitHub repo, as the reviewer asked.

## Already solid (no change needed)

- RPD definition uses BKS across all methods (exact + heuristic) — matches Q2.
- Nonparametric pipeline (Shapiro–Wilk → Friedman → Bonferroni Wilcoxon + effect
  sizes) — correct structure for Q3.
- Per-instance avg/std reported, not just best-so-far — core of Q2.

## Decisions needed from author before editing

1. How many runs were actually performed per method? (drives the whole Q3 answer)
2. Was SSM-SA single-run or averaged? (drives critical item #2)
3. Is greedy-baseline RPD data available to add as a column?
