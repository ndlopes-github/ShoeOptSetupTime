# Supplementary Material

Reproducibility artifacts for the paper *Parallel machine scheduling with
simultaneous interruptions: a case study in shoe sole production*
(J. Orestes Cerdeira, Ricardo Enguiça, Nuno Lopes).

This directory serves **two purposes**:

1. **Statistical analysis requested by a reviewer** — the per-instance
   *average* RPD (Relative Percentage Deviation) of every stochastic method
   relative to the best-known solution (taken over all methods, exact and
   heuristic), together with a formal nonparametric analysis (Shapiro–Wilk
   normality, Friedman omnibus test, and Bonferroni-corrected pairwise Wilcoxon
   signed-rank tests with effect sizes), following Molina et al. (2021).
2. **irace parameter-tuning scripts** — the configurations used to calibrate
   each heuristic. These also fix the number of independent runs (`N_runs`) per
   method, which is exactly the number of runs averaged in the statistical
   analysis (`N_runs` = 181/196 for SA at β = 3/6, 144 for GRASP, 100 for GA;
   SSM-SA is single-run by design, as each run entails repeated MILP solves).

## Directory layout

| Path | Contents |
|------|----------|
| `GA/`, `GRASP/`, `MILP/`, `SA/`, `SSM-SA/` | Per-method result archives (`.tar.bz2`): per-run logs (and, for MILP and SSM-SA, solution matrices) for β = 3 and β = 6. The per-run logs are the raw inputs for the average-RPD computation. |
| `GREEDY/` | Reserved for greedy-heuristic outputs; the greedy method is deterministic (single solution per instance) and is reported directly in the paper. |
| `INSTANCES/MILP/` | `milp.tar.bz2` — instances for the exact solver. |
| `INSTANCES/HEURISTICS/` | `heuristics.tar.bz2` — instances for the heuristic methods. |
| `IRACE_SCRIPTS/` | irace tuning setup, one directory per method (`irace_ga/`, `irace_grasp/`, `irace_sa/`, `irace_ssm_sa/`): `scenario.txt`, `parameters.txt`, `forbidden.txt`, `target-runner`, the instance list, and the per-β tuning logs (`logs_beta3/`, `logs_beta6/`). The selected configurations — including `N_runs` — are reported in the *Computational experiments* section of the paper. |
| `STATS_SCRIPTS/` | Statistical analysis tools and their outputs (see below). |
| [`rpd_boxplot.png`](rpd_boxplot.png) | Boxplot of the per-instance average RPD distributions across the four heuristics. |
| [`supplementary_rpd.md`](supplementary_rpd.md) | **Per-instance average RPD table for all methods and the full statistical report** — the main analysis document referenced from the paper. |

## `STATS_SCRIPTS/`

Julia scripts (run with `julia <script>.jl` from inside `STATS_SCRIPTS/`) and
their committed outputs.

| File | Role |
|------|------|
| `compute_rpd.jl` | Computes the best-known solution (BKS) per instance over all methods, the per-run RPD for SA/GRASP/GA and the single-run RPD for SSM-SA, and writes `rpd_results.csv` and the per-instance table. Requires the per-method run archives to be extracted (see the header of the script for the `RUN_DIR` setting). |
| `statistical_analysis.jl` | Shapiro–Wilk normality, Friedman omnibus test, and Bonferroni-corrected pairwise Wilcoxon signed-rank tests with rank-biserial effect sizes; produces the boxplot. Reads `rpd_results.csv`. |
| `reviewer_stats.jl` | Run-count summary, matched-set (109-instance) RPD summary, and matched-set Shapiro–Wilk normality consistent with the Friedman/Wilcoxon instance set. Reads `rpd_results.csv`. |
| `rpd_results.csv` | Per-instance RPD table (input to the analysis scripts). |
| `tables_04_to_13_combined.csv` | MILP and SSM-SA per-instance objective values used to build `rpd_results.csv`. |
| `stats_normality.csv`, `stats_friedman.txt`, `stats_wilcoxon.csv` | Outputs of `statistical_analysis.jl`. |
| `stats_runcounts.csv`, `stats_summary_full.csv`, `stats_summary_matched.csv`, `stats_normality_matched.csv` | Outputs of `reviewer_stats.jl`. |
| `rpd_boxplot.pdf`, `rpd_boxplot.png` | Boxplot figure (PDF/PNG). |

### Required Julia packages

```julia
import Pkg
Pkg.add(["CSV", "DataFrames", "HypothesisTests", "StatsBase",
         "Distributions", "Statistics", "StatsPlots", "Plots"])
```

### Reproducing the statistical analysis

```bash
cd STATS_SCRIPTS
# (1) regenerate rpd_results.csv from the run archives — see compute_rpd.jl header
julia compute_rpd.jl
# (2) omnibus + post-hoc tests and boxplot
julia statistical_analysis.jl
# (3) run counts, matched-set summary, and matched-set normality
julia reviewer_stats.jl
```

Steps (2) and (3) run directly from the committed `rpd_results.csv`; step (1)
is only needed to rebuild that table from the raw per-run logs.
