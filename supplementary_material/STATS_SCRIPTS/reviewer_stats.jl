#!/usr/bin/env julia
"""
reviewer_stats.jl — Supplementary statistics requested during peer review.

Complements compute_rpd.jl and statistical_analysis.jl. It reads the per-instance
RPD table (rpd_results.csv) and produces three things that the existing scripts
do not, all needed to answer reviewer comments #2 and #3:

  1. Run-count summary: the number of independent runs averaged per method
     (N_runs), per value of beta — confirming the >= 30 independent runs
     requested by the reviewer (SSM-SA is single-run by design; see paper).
  2. Matched-set summary: mean/std/median/min/max RPD computed BOTH over the
     full set available per method AND over the 109-instance set on which all
     four heuristics are applicable, so the headline means are comparable.
  3. Matched-set Shapiro-Wilk normality (n = 109), consistent with the
     instance set actually used by the Friedman and Wilcoxon tests.

Required packages (install once):
    import Pkg; Pkg.add(["CSV", "DataFrames", "HypothesisTests", "Statistics"])

Usage:
    julia reviewer_stats.jl

Outputs (written next to this script):
    stats_runcounts.csv         — N_runs per method per beta
    stats_summary_full.csv      — summary RPD over each method's full set
    stats_summary_matched.csv   — summary RPD over the common 109-instance set
    stats_normality_matched.csv — Shapiro-Wilk on the 109 matched instances
"""

using CSV
using DataFrames
using HypothesisTests
using Statistics
using Printf

const SCRIPT_DIR = dirname(abspath(@__FILE__))
const IN_CSV = joinpath(SCRIPT_DIR, "rpd_results.csv")
const OUT_RUNS = joinpath(SCRIPT_DIR, "stats_runcounts.csv")
const OUT_FULL = joinpath(SCRIPT_DIR, "stats_summary_full.csv")
const OUT_MATCHED = joinpath(SCRIPT_DIR, "stats_summary_matched.csv")
const OUT_NORM_MATCHED = joinpath(SCRIPT_DIR, "stats_normality_matched.csv")

# Methods and the rpd_results.csv column holding their per-instance (avg) RPD.
const METHODS = [
    ("SA", :SA_avg_RPD),
    ("GRASP", :GRASP_avg_RPD),
    ("GA", :GA_avg_RPD),
    ("SSM-SA", :SSM_SA_RPD),
    ("MILP", :MILP_RPD),
]

cleanvals(col) = Float64[x for x in col if !ismissing(x)]

function summarise(vals::Vector{Float64})
    isempty(vals) && return nothing
    (n=length(vals), mean=mean(vals), std=(length(vals) < 2 ? NaN : std(vals)),
     median=median(vals), min=minimum(vals), max=maximum(vals))
end

function summary_frame(df::DataFrame)
    rows = NamedTuple[]
    for (label, col) in METHODS
        s = summarise(cleanvals(df[!, col]))
        s === nothing && continue
        push!(rows, (Algorithm=label, n=s.n, mean_RPD=s.mean, std_RPD=s.std,
                     median_RPD=s.median, min_RPD=s.min, max_RPD=s.max))
    end
    DataFrame(rows)
end

function print_summary(title, sdf)
    println("\n=== $title ===")
    @printf("%-8s %4s %9s %9s %9s %9s %9s\n",
            "Algo", "n", "mean", "std", "median", "min", "max")
    println(repeat('-', 62))
    for r in eachrow(sdf)
        @printf("%-8s %4d %9.4f %9.4f %9.4f %9.4f %9.4f\n",
                r.Algorithm, r.n, r.mean_RPD, r.std_RPD, r.median_RPD, r.min_RPD, r.max_RPD)
    end
end

function main()
    df = CSV.read(IN_CSV, DataFrame; missingstring=["", "—"])
    println("Loaded $(nrow(df)) instances from $(basename(IN_CSV)).")

    # ── 1. Run-count summary (per method, per beta) ────────────────────────────
    runcols = [("SA", :SA_n), ("GRASP", :GRASP_n), ("GA", :GA_n)]
    runrows = NamedTuple[]
    for b in sort(unique(df.beta)), (label, col) in runcols
        vals = df[df.beta .== b, col]
        push!(runrows, (Algorithm=label, beta=b,
                        N_runs_min=minimum(vals), N_runs_max=maximum(vals)))
    end
    # SSM-SA is single-run by design.
    for b in sort(unique(df.beta))
        push!(runrows, (Algorithm="SSM-SA", beta=b, N_runs_min=1, N_runs_max=1))
    end
    runs_df = DataFrame(runrows)
    CSV.write(OUT_RUNS, runs_df)
    println("\n=== Independent runs averaged per instance (N_runs) ===")
    for r in eachrow(runs_df)
        rng = r.N_runs_min == r.N_runs_max ? "$(r.N_runs_min)" :
              "$(r.N_runs_min)-$(r.N_runs_max)"
        @printf("  %-8s  beta=%d  N_runs=%s\n", r.Algorithm, r.beta, rng)
    end
    println("Written: $OUT_RUNS")

    # ── 2. Summaries: full set vs matched 109-instance set ─────────────────────
    full_sdf = summary_frame(df)
    matched = df[.!ismissing.(df.SSM_SA_RPD), :]
    matched_sdf = summary_frame(matched)

    CSV.write(OUT_FULL, full_sdf)
    CSV.write(OUT_MATCHED, matched_sdf)
    print_summary("Summary RPD — full set (per-method availability)", full_sdf)
    print_summary("Summary RPD — matched set (n=$(nrow(matched)) common instances)", matched_sdf)
    println("\nWritten: $OUT_FULL")
    println("Written: $OUT_MATCHED")

    # ── 3. Matched-set Shapiro-Wilk normality (n=109) ──────────────────────────
    println("\n=== Shapiro-Wilk normality on the $(nrow(matched)) matched instances ===")
    normrows = NamedTuple[]
    for (label, col) in METHODS[1:4]   # four heuristics only
        v = cleanvals(matched[!, col])
        sw = ShapiroWilkTest(v)
        p = pvalue(sw)
        @printf("  %-8s n=%3d  W=%.4f  p=%.4e  %s\n",
                label, length(v), sw.W, p, p ≥ 0.05 ? "normal" : "NON-NORMAL")
        push!(normrows, (Algorithm=label, n=length(v), W=sw.W, p_value=p,
                         Normal_p_ge_0_05=(p ≥ 0.05)))
    end
    CSV.write(OUT_NORM_MATCHED, DataFrame(normrows))
    println("Written: $OUT_NORM_MATCHED")
end

main()
