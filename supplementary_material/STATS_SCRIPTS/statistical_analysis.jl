#!/usr/bin/env julia
"""
statistical_analysis.jl — Statistical comparison of heuristics via per-instance average RPD.

Julia port of statistical_analysis.py.

Steps (following Derrac et al., Swarm Evol. Comput. 1(1), 3–18, 2011):
  1. Shapiro–Wilk normality test per algorithm
  2. Friedman test (non-parametric overall comparison)
  3. Wilcoxon signed-rank pairwise post-hoc tests + Bonferroni correction
  4. Rank-biserial correlation (effect size) for each pair
  5. Boxplot of per-instance average RPD distributions

NOTE: SSM-SA has only 1 run per instance; its RPD is a single observation
(no within-instance variance). Acknowledged in the statistical analysis.
MILP is excluded: it always achieves RPD=0 on proven-optimal instances and
is not defined on others.

Required packages (install once):
    import Pkg
    Pkg.add(["CSV", "DataFrames", "HypothesisTests", "StatsBase",
             "Distributions", "StatsPlots", "Plots"])

Usage:
    julia statistical_analysis.jl

Outputs:
    <script_dir>/stats_normality.csv   — Shapiro–Wilk results
    <script_dir>/stats_friedman.txt    — Friedman test result
    <script_dir>/stats_wilcoxon.csv    — Pairwise Wilcoxon + effect sizes
    <script_dir>/rpd_boxplot.pdf       — Boxplot figure (PDF)
    <script_dir>/rpd_boxplot.png       — Boxplot figure (PNG)
"""

using CSV
using DataFrames
using HypothesisTests
using StatsBase: tiedrank
using Distributions: Chisq, ccdf
using Statistics
using Printf
using Plots
using StatsPlots

# ── Paths ──────────────────────────────────────────────────────────────────────
const SCRIPT_DIR = dirname(abspath(@__FILE__))
const IN_CSV = joinpath(SCRIPT_DIR, "rpd_results.csv")
const OUT_NORM = joinpath(SCRIPT_DIR, "stats_normality.csv")
const OUT_FRIED = joinpath(SCRIPT_DIR, "stats_friedman.txt")
const OUT_WILCOX = joinpath(SCRIPT_DIR, "stats_wilcoxon.csv")
const OUT_PLOT = joinpath(SCRIPT_DIR, "rpd_boxplot.pdf")
const OUT_PLOT_PNG = joinpath(SCRIPT_DIR, "rpd_boxplot.png")

# ── Load data ──────────────────────────────────────────────────────────────────
"""
Returns:
  data          — Dict{String, Vector{Float64}}: all per-instance values per algo
                  (instances missing for that algo are excluded)
  complete_rows — Vector{Dict}: rows where ALL four heuristics have a value
                  (used for Friedman and Wilcoxon)
"""
function load_rpd(path::String)
    algo_cols = ["SA_avg_RPD", "GRASP_avg_RPD", "GA_avg_RPD", "SSM_SA_RPD"]
    df = CSV.read(path, DataFrame; missingstring=["", "—"])

    data = Dict(a => Float64[] for a in algo_cols)
    complete_rows = Vector{Dict{String,Float64}}()

    for row in eachrow(df)
        vals = Dict{String,Union{Float64,Missing}}()
        for a in algo_cols
            vals[a] = ismissing(row[Symbol(a)]) ? missing : Float64(row[Symbol(a)])
        end
        for a in algo_cols
            ismissing(vals[a]) || push!(data[a], vals[a])
        end
        if !any(ismissing, values(vals))
            push!(complete_rows, Dict(a => vals[a] for a in algo_cols))
        end
    end

    return data, complete_rows
end

# ── Helpers ────────────────────────────────────────────────────────────────────
"""
Rank-biserial correlation for the Wilcoxon signed-rank test.
Formula: r = 1 − (2W) / (n(n+1)/2)
where W = min(W⁺, W⁻) over non-zero pairwise differences.
"""
function rank_biserial(x::Vector{Float64}, y::Vector{Float64})::Float64
    diffs = [xi - yi for (xi, yi) in zip(x, y) if xi ≠ yi]
    n = length(diffs)
    n == 0 && return 0.0

    abs_d = abs.(diffs)
    ranks = tiedrank(abs_d)          # handles ties; ranks start at 1
    W_plus = sum(ranks[i] for i in 1:n if diffs[i] > 0; init=0.0)
    W_minus = sum(ranks[i] for i in 1:n if diffs[i] < 0; init=0.0)
    W = min(W_plus, W_minus)
    max_W = n * (n + 1) / 2
    return 1.0 - (2.0 * W) / max_W
end

function effect_magnitude(r::Float64)::String
    ar = abs(r)
    ar < 0.1 && return "negligible"
    ar < 0.3 && return "small"
    ar < 0.5 && return "medium"
    return "large"
end

fmt(::Missing, digits::Int=4) = "—"
fmt(v::Real, digits::Int=4) = @sprintf("%.*f", digits, Float64(v))

# ── Friedman test (manual implementation) ─────────────────────────────────────
"""
Non-parametric Friedman test on k repeated-measures groups.

Arguments:
  groups — k vectors of length n (same n for all)

Returns (χ², p-value) where χ² ~ Chisq(k-1) under H₀.
"""
function friedman_test(groups::Vector{Vector{Float64}})
    k = length(groups)
    n = length(groups[1])
    @assert all(length(g) == n for g in groups) "All groups must have equal length."

    # Rank each block (row) across treatments
    # ranked[i, j] = rank of treatment j in block i
    ranked = Matrix{Float64}(undef, n, k)
    for i in 1:n
        row = Float64[groups[j][i] for j in 1:k]
        ranked[i, :] = tiedrank(row)
    end

    # Column rank sums
    R = vec(sum(ranked; dims=1))

    # Friedman statistic Q ~ χ²(k-1) for large n
    Q = 12.0 / (n * k * (k + 1)) * sum(R .^ 2) - 3.0 * n * (k + 1)
    p = ccdf(Chisq(k - 1), Q)
    return Q, p
end

# ── Boxplot ─────────────────────────────────────────────────────────────────────
function plot_boxplot(data::Dict{String,Vector{Float64}}, algo_cols, algo_labels)
    try

        gr()
        plot_data = [data[c] for c in algo_cols]
        colors = ["#4C72B0" "#DD8452" "#55A868" "#C44E52"]
        p = StatsPlots.boxplot(
            reshape(collect(algo_labels), 1, :), plot_data;
            fillalpha=0.7,
            color=colors,
            linecolor=:black,
            whisker_width=0.5,
            outliers=true,
            ylabel="Average RPD (%)",
            xlabel="Algorithm",
            title="Distribution of per-instance average RPD",
            legend=false,
            size=(640, 420),
        )
        Plots.savefig(p, OUT_PLOT)
        println("Written: $OUT_PLOT")
        Plots.savefig(p, OUT_PLOT_PNG)
        println("Written: $OUT_PLOT_PNG")
    catch e
        @warn "Could not produce boxplot (StatsPlots/Plots not available): $e"
    end
end

# ── Main ───────────────────────────────────────────────────────────────────────
function main()
    data, complete_rows = load_rpd(IN_CSV)

    algo_cols = ["SA_avg_RPD", "GRASP_avg_RPD", "GA_avg_RPD", "SSM_SA_RPD"]
    algo_labels = ["SA", "GRASP", "GA", "SSM-SA"]

    println("Loaded data: $(map(a -> length(data[a]), algo_cols)) values per algo")
    println("Complete rows (all 4 algos present): $(length(complete_rows))\n")

    # ── 1. Shapiro–Wilk normality test ────────────────────────────────────────
    println("=== Shapiro–Wilk Normality Test ===")
    norm_rows = NamedTuple[]
    for (col, label) in zip(algo_cols, algo_labels)
        v = data[col]
        sw = ShapiroWilkTest(v)
        W_stat = sw.W
        p_val = pvalue(sw)
        normal = p_val ≥ 0.05
        @printf("  %-8s  n=%3d  W=%.4f  p=%.4e  %s\n",
            label, length(v), W_stat, p_val, normal ? "normal" : "NON-NORMAL")
        push!(norm_rows, (
            Algorithm=label,
            n=length(v),
            W=W_stat,
            p_value=p_val,
            Normal_p_ge_0_05=normal,
        ))
    end
    CSV.write(OUT_NORM, DataFrame(norm_rows))
    println("\nWritten: $OUT_NORM")

    # ── 2. Friedman test ──────────────────────────────────────────────────────
    n_complete = length(complete_rows)
    groups = [Float64[r[c] for r in complete_rows] for c in algo_cols]
    Q, fp = friedman_test(groups)
    sig = fp < 0.05
    println("\n=== Friedman Test (n=$n_complete complete instances) ===")
    @printf("  χ² = %.4f,  p = %.4e\n", Q, fp)
    println("  $(sig ? "Significant difference detected." : "No significant difference.")")

    open(OUT_FRIED, "w") do io
        println(io, "Friedman test on $n_complete instances with all 4 heuristics")
        @printf(io, "chi2 = %.6f\n", Q)
        @printf(io, "p    = %.6e\n", fp)
        println(io, "Significant (p<0.05): $sig")
    end
    println("Written: $OUT_FRIED")

    # ── 3. Wilcoxon signed-rank pairwise post-hoc ─────────────────────────────
    pairs = [(i, j) for i in 1:length(algo_cols) for j in (i+1):length(algo_cols)]
    n_pairs = length(pairs)
    α_bonf = 0.05 / n_pairs
    println("\n=== Wilcoxon Signed-Rank (pairwise, Bonferroni α=$(round(α_bonf; digits=4))) ===")

    wilcox_rows = NamedTuple[]
    for (i, j) in pairs
        xi = Float64[r[algo_cols[i]] for r in complete_rows]
        xj = Float64[r[algo_cols[j]] for r in complete_rows]

        # HypothesisTests for two-sided p-value
        st = SignedRankTest(xi, xj)
        p_raw = pvalue(st; tail=:both)
        p_bonf = min(p_raw * n_pairs, 1.0)

        # W statistic: min(T⁺, T⁻), consistent with scipy.stats.wilcoxon
        T_plus = st.W
        n_nonzero = count(!iszero, xi .- xj)
        T_total = n_nonzero * (n_nonzero + 1) / 2
        W_stat = min(T_plus, T_total - T_plus)

        # Effect size
        r_rb = rank_biserial(xi, xj)
        mag = effect_magnitude(r_rb)

        sig_raw = p_raw < 0.05
        sig_bonf = p_bonf < 0.05

        @printf("  %-8s vs %-8s  W=%7.1f  p=%10.4e  p_bonf=%10.4e  r=%+.3f (%s)  %s %s\n",
            algo_labels[i], algo_labels[j],
            W_stat, p_raw, p_bonf,
            r_rb, mag,
            sig_raw ? "*" : " ",
            sig_bonf ? "*" : " ")

        push!(wilcox_rows, (
            Algo_A=algo_labels[i],
            Algo_B=algo_labels[j],
            W_stat=W_stat,
            p_value=p_raw,
            p_bonferroni=p_bonf,
            sig_raw=sig_raw,
            sig_bonferroni=sig_bonf,
            rank_biserial=r_rb,
            effect_size=mag,
        ))
    end
    CSV.write(OUT_WILCOX, DataFrame(wilcox_rows))
    println("\nWritten: $OUT_WILCOX")

    # ── 4. Boxplot ────────────────────────────────────────────────────────────
    plot_boxplot(data, algo_cols, algo_labels)
end

main()
