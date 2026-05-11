#!/usr/bin/env python3
"""
statistical_analysis.py — Statistical comparison of heuristics via per-instance average RPD.

Steps (following Derrac et al., Swarm Evol. Comput. 64, 100888):
  1. Shapiro-Wilk normality test per algorithm
  2. Friedman test (non-parametric overall comparison)
  3. Wilcoxon signed-rank pairwise post-hoc tests + Bonferroni correction
  4. Rank-biserial correlation (effect size) for each pair
  5. Boxplot of per-instance average RPD distributions

NOTE: SSM-SA has only 1 run per instance; its RPD is a single observation per instance
(no within-instance variance). This is valid for a between-instance comparison but is
acknowledged in a footnote.
MILP is excluded from the statistical comparison: it always achieves RPD=0 on proven-optimal
instances and is not defined on non-optimal ones — it is not a heuristic.

Outputs:
  aux_temp/scripts/stats_normality.csv   — Shapiro-Wilk results
  aux_temp/scripts/stats_friedman.txt    — Friedman test result
  aux_temp/scripts/stats_wilcoxon.csv    — Pairwise Wilcoxon + effect sizes
  aux_temp/scripts/rpd_boxplot.pdf       — Boxplot figure
"""

import csv
import itertools
import math
from pathlib import Path

import scipy.stats as stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
IN_CSV     = SCRIPT_DIR / "rpd_results.csv"
OUT_NORM   = SCRIPT_DIR / "stats_normality.csv"
OUT_FRIED  = SCRIPT_DIR / "stats_friedman.txt"
OUT_WILCOX = SCRIPT_DIR / "stats_wilcoxon.csv"
OUT_PLOT   = SCRIPT_DIR / "rpd_boxplot.pdf"

# ── Load data ─────────────────────────────────────────────────────────────────
def load_rpd(path):
    """
    Returns dict: algo -> list of per-instance average RPD values (float).
    Instances where a value is missing (empty string) are excluded per-algo.
    Only instances where ALL four heuristics have data are used for Friedman/Wilcoxon.
    """
    algos = ['SA_avg_RPD', 'GRASP_avg_RPD', 'GA_avg_RPD', 'SSM_SA_RPD']
    data  = {a: [] for a in algos}
    complete_rows = []   # rows where all 4 algos have data

    with open(path, newline='') as f:
        for row in csv.DictReader(f):
            vals = {}
            for a in algos:
                v = row[a].strip()
                vals[a] = float(v) if v and v != '—' else None
            for a in algos:
                if vals[a] is not None:
                    data[a].append(vals[a])
            if all(v is not None for v in vals.values()):
                complete_rows.append({a: vals[a] for a in algos})

    return data, complete_rows


# ── Helpers ───────────────────────────────────────────────────────────────────
def rank_biserial(x, y):
    """
    Rank-biserial correlation for Wilcoxon signed-rank test.
    r = 1 - (2 * W) / (n * (n+1)/2)  where W is the Wilcoxon statistic.
    """
    diffs = [xi - yi for xi, yi in zip(x, y) if (xi - yi) != 0]
    n = len(diffs)
    if n == 0:
        return 0.0
    abs_ranks = sorted(range(n), key=lambda i: abs(diffs[i]))
    # rank starting at 1
    ranks = [0] * n
    i = 0
    while i < n:
        j = i
        while j < n - 1 and abs(diffs[abs_ranks[j + 1]]) == abs(diffs[abs_ranks[i]]):
            j += 1
        avg_rank = (i + 1 + j + 1) / 2
        for k in range(i, j + 1):
            ranks[abs_ranks[k]] = avg_rank
        i = j + 1

    W_plus  = sum(ranks[i] for i in range(n) if diffs[i] > 0)
    W_minus = sum(ranks[i] for i in range(n) if diffs[i] < 0)
    W = min(W_plus, W_minus)
    max_W = n * (n + 1) / 2
    r = 1 - (2 * W) / max_W
    return r


def effect_magnitude(r):
    ar = abs(r)
    if ar < 0.1:   return "negligible"
    if ar < 0.3:   return "small"
    if ar < 0.5:   return "medium"
    return "large"


def fmt(v, digits=4):
    return '—' if v is None else f'{v:.{digits}f}'


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    data, complete_rows = load_rpd(IN_CSV)

    algo_cols   = ['SA_avg_RPD', 'GRASP_avg_RPD', 'GA_avg_RPD', 'SSM_SA_RPD']
    algo_labels = ['SA',         'GRASP',          'GA',          'SSM-SA']

    print(f"Loaded data: {[len(data[a]) for a in algo_cols]} values per algo")
    print(f"Complete rows (all 4 algos): {len(complete_rows)}\n")

    # ── 1. Shapiro-Wilk ────────────────────────────────────────────────────────
    print("=== Shapiro-Wilk Normality Test ===")
    norm_rows = []
    for col, label in zip(algo_cols, algo_labels):
        v = data[col]
        stat, p = stats.shapiro(v)
        normal = p >= 0.05
        print(f"  {label:8s}  n={len(v):3d}  W={stat:.4f}  p={p:.4e}  {'normal' if normal else 'NON-NORMAL'}")
        norm_rows.append({'Algorithm': label, 'n': len(v), 'W': stat, 'p_value': p,
                          'Normal (p>=0.05)': normal})

    with open(OUT_NORM, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(norm_rows[0].keys()))
        w.writeheader()
        w.writerows(norm_rows)
    print(f"\nWritten: {OUT_NORM}")

    # ── 2. Friedman test (on complete rows only) ───────────────────────────────
    groups = [[r[c] for r in complete_rows] for c in algo_cols]
    fstat, fp = stats.friedmanchisquare(*groups)
    print(f"\n=== Friedman Test (n={len(complete_rows)} complete instances) ===")
    print(f"  χ² = {fstat:.4f},  p = {fp:.4e}")
    significant = fp < 0.05
    print(f"  {'Significant difference detected.' if significant else 'No significant difference.'}")

    with open(OUT_FRIED, 'w') as f:
        f.write(f"Friedman test on {len(complete_rows)} instances with all 4 heuristics\n")
        f.write(f"chi2 = {fstat:.6f}\n")
        f.write(f"p    = {fp:.6e}\n")
        f.write(f"Significant (p<0.05): {significant}\n")
    print(f"Written: {OUT_FRIED}")

    # ── 3. Wilcoxon pairwise post-hoc ─────────────────────────────────────────
    pairs = list(itertools.combinations(range(len(algo_cols)), 2))
    n_pairs = len(pairs)
    print(f"\n=== Wilcoxon Signed-Rank (pairwise, Bonferroni α={0.05/n_pairs:.4f}) ===")

    wilcox_rows = []
    for i, j in pairs:
        xi = [r[algo_cols[i]] for r in complete_rows]
        xj = [r[algo_cols[j]] for r in complete_rows]
        wstat, wp = stats.wilcoxon(xi, xj, alternative='two-sided')
        wp_bonf = min(wp * n_pairs, 1.0)
        r_rb    = rank_biserial(xi, xj)
        mag     = effect_magnitude(r_rb)
        sig_raw  = wp      < 0.05
        sig_bonf = wp_bonf < 0.05
        print(f"  {algo_labels[i]:8s} vs {algo_labels[j]:8s}  "
              f"W={wstat:.1f}  p={wp:.4e}  p_bonf={wp_bonf:.4e}  "
              f"r={r_rb:+.3f} ({mag})  "
              f"{'*' if sig_raw else ' '} {'*' if sig_bonf else ' '}")
        wilcox_rows.append({
            'Algo_A':        algo_labels[i],
            'Algo_B':        algo_labels[j],
            'W_stat':        wstat,
            'p_value':       wp,
            'p_bonferroni':  wp_bonf,
            'sig_raw':       sig_raw,
            'sig_bonferroni': sig_bonf,
            'rank_biserial': r_rb,
            'effect_size':   mag,
        })

    with open(OUT_WILCOX, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(wilcox_rows[0].keys()))
        w.writeheader()
        w.writerows(wilcox_rows)
    print(f"\nWritten: {OUT_WILCOX}")

    # ── 4. Boxplot ─────────────────────────────────────────────────────────────
    plot_data   = [data[c]   for c in algo_cols]
    plot_labels = algo_labels

    fig, ax = plt.subplots(figsize=(7, 4.5))
    bp = ax.boxplot(plot_data, labels=plot_labels, patch_artist=True,
                    medianprops=dict(color='black', linewidth=1.5),
                    flierprops=dict(marker='o', markersize=3, alpha=0.5))

    colors = ['#4C72B0', '#DD8452', '#55A868', '#C44E52']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_ylabel('Average RPD (%)', fontsize=11)
    ax.set_xlabel('Algorithm', fontsize=11)
    ax.set_title('Distribution of per-instance average RPD', fontsize=11)
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.grid(axis='y', linestyle='--', alpha=0.4)
    plt.tight_layout()
    fig.savefig(OUT_PLOT, dpi=150)
    print(f"Written: {OUT_PLOT}")


if __name__ == '__main__':
    main()
