#!/usr/bin/env python3
"""
compute_rpd.py — Compute per-run RPD for SA, GRASP, GA and single-value RPD for SSM-SA.

MILP and SSM-SA costs are read from the paper tables in main.tex:
  - Plain number          → optimal / single-run cost
  - "B xxx (F yyy)"      → MILP did not find optimal; excluded from BKS and RPD (marked —)
  - "**"                 → SSM-SA not run for this instance (marked —)

BKS per instance = min over all available feasible costs.

Outputs:
  aux_temp/scripts/rpd_results.csv   — per-instance RPD stats (used by statistical_analysis.py)
  supplementary_rpd.md               — GitHub supplementary material (repo root)
"""

import re
import csv
import math
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT     = Path(__file__).resolve().parent.parent.parent   # repo root
AUX      = ROOT / "aux_temp"
MAIN_TEX = ROOT / "main.tex"
OUT_CSV  = Path(__file__).resolve().parent / "rpd_results.csv"
OUT_MD   = ROOT / "supplementary_rpd.md"

# ── Regex ─────────────────────────────────────────────────────────────────────
RE_ORDER   = re.compile(r'\\multirow\{4\}\{\*\}\{\$O_(\d)\$\}')
RE_BF      = re.compile(r'B\s*[\d.]+\s*\(F\s*[\d.]+\)')
RE_STRUCT  = re.compile(
    r'\\(?:toprule|bottomrule|midrule|cmidrule\s*\([^)]*\)\s*\{[^}]*\})'
)

# ── Cell parsing ──────────────────────────────────────────────────────────────
def parse_t_cell(raw):
    """
    Parse a t-column cell.
    Returns float, or None if B/F (non-optimal MILP) or missing (**).
    """
    c = raw.strip()
    if not c or c in ('**', '*', '--'):
        return None
    if RE_BF.search(c):
        return None          # non-optimal MILP — exclude from BKS
    try:
        return float(c)
    except ValueError:
        return None

# ── Parse LaTeX tables ────────────────────────────────────────────────────────
def parse_latex_tables(tex_path):
    """
    Returns dict keyed by (scenario: int, beta: int, order: int, p: int)
    with values {'milp': float|None, 'ssm_sa': float|None}.
    """
    tex = Path(tex_path).read_text(encoding='utf-8')
    # Split on table environments
    blocks = re.split(r'(?=\\begin\{table\})', tex)

    result = {}
    for block in blocks:
        label_m = re.search(
            r'\\label\{table:combined_methods(_beta6)?_(\d)\}', block
        )
        if not label_m:
            continue
        beta     = 6 if label_m.group(1) else 3
        scenario = int(label_m.group(2))

        # Extract tabular content
        tab_m = re.search(
            r'\\begin\{tabular\}.*?\\end\{tabular\}', block, re.DOTALL
        )
        if not tab_m:
            continue

        tabular = tab_m.group(0)
        # Split rows on \\
        rows = re.split(r'\\\\', tabular)

        current_order = None
        for row in rows:
            # Detect order from \multirow{4}
            order_m = RE_ORDER.search(row)
            if order_m:
                current_order = int(order_m.group(1))

            # Strip structural commands, then check for cell data
            cleaned = RE_STRUCT.sub('', row)
            cells = cleaned.split('&')
            if len(cells) < 4:
                continue

            # cells[1] = p value
            p_str = cells[1].strip()
            try:
                p = int(p_str)
            except ValueError:
                continue
            if p not in (2, 3, 4, 5):
                continue

            if current_order is None:
                continue

            milp_t   = parse_t_cell(cells[2])
            ssm_sa_t = parse_t_cell(cells[3])

            key = (scenario, beta, current_order, p)
            result[key] = {'milp': milp_t, 'ssm_sa': ssm_sa_t}

    return result

# ── Load heuristic run CSVs ───────────────────────────────────────────────────
def load_csv_costs(algo, beta, scenario, order, p):
    """Return list of float costs from the appropriate run CSV."""
    folder_map = {
        'sa':    AUX / f'logs_sa_beta{beta}'    / 'sa_runs',
        'grasp': AUX / f'logs_grasp_beta{beta}' / 'grasp_runs',
        'ga':    AUX / f'logs_ga_beta{beta}'    / 'ga_runs',
    }
    fpath = folder_map[algo] / f'H_O{order}_#{scenario}_{p}p.csv'
    if not fpath.exists():
        print(f"  WARNING: missing CSV {fpath}")
        return []
    costs = []
    with open(fpath, newline='') as f:
        for row in csv.DictReader(f):
            try:
                costs.append(float(row['cost']))
            except (KeyError, ValueError):
                pass
    return costs

# ── Statistics helpers ────────────────────────────────────────────────────────
def _vals(seq):
    return [v for v in seq if v is not None]

def mean(seq):
    v = _vals(seq)
    return sum(v) / len(v) if v else None

def std(seq):
    v = _vals(seq)
    if len(v) < 2:
        return None
    m = mean(v)
    return math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1))

def median(seq):
    v = sorted(_vals(seq))
    n = len(v)
    if not n:
        return None
    mid = n // 2
    return v[mid] if n % 2 else (v[mid - 1] + v[mid]) / 2

def rpd(cost, bks):
    if cost is None or bks is None or bks == 0:
        return None
    return (cost - bks) / bks * 100.0

def fmt(v, digits=4):
    return '—' if v is None else f'{v:.{digits}f}'

# ── Summary helper ────────────────────────────────────────────────────────────
def summarise(records, col):
    v = _vals([r[col] for r in records])
    if not v:
        return None
    return {
        'n':      len(v),
        'mean':   mean(v),
        'std':    std(v),
        'median': median(v),
        'min':    min(v),
        'max':    max(v),
    }

# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    print("Parsing LaTeX tables …")
    table_values = parse_latex_tables(MAIN_TEX)
    print(f"  Found {len(table_values)} instance entries in tables.")

    records = []
    missing_table = []

    for beta in (3, 6):
        for scenario in range(1, 6):
            for order in (1, 2, 3):
                for p in (2, 3, 4, 5):
                    key = (scenario, beta, order, p)
                    tv  = table_values.get(key)
                    if tv is None:
                        missing_table.append(key)
                        tv = {'milp': None, 'ssm_sa': None}

                    milp_t   = tv['milp']
                    ssm_sa_t = tv['ssm_sa']

                    sa_costs    = load_csv_costs('sa',    beta, scenario, order, p)
                    grasp_costs = load_csv_costs('grasp', beta, scenario, order, p)
                    ga_costs    = load_csv_costs('ga',    beta, scenario, order, p)

                    # BKS: min of all available feasible costs
                    all_costs = []
                    if milp_t   is not None: all_costs.append(milp_t)
                    if ssm_sa_t is not None: all_costs.append(ssm_sa_t)
                    all_costs.extend(sa_costs)
                    all_costs.extend(grasp_costs)
                    all_costs.extend(ga_costs)

                    if not all_costs:
                        print(f"  WARNING: no cost data for {key}, skipping.")
                        continue

                    bks = min(all_costs)

                    sa_rpds    = [rpd(c, bks) for c in sa_costs]
                    grasp_rpds = [rpd(c, bks) for c in grasp_costs]
                    ga_rpds    = [rpd(c, bks) for c in ga_costs]

                    records.append({
                        'instance':      f'H_O{order}_#{scenario}_{p}p',
                        'beta':          beta,
                        'BKS':           bks,
                        'SA_avg_RPD':    mean(sa_rpds),
                        'SA_std_RPD':    std(sa_rpds),
                        'SA_n':          len(sa_costs),
                        'GRASP_avg_RPD': mean(grasp_rpds),
                        'GRASP_std_RPD': std(grasp_rpds),
                        'GRASP_n':       len(grasp_costs),
                        'GA_avg_RPD':    mean(ga_rpds),
                        'GA_std_RPD':    std(ga_rpds),
                        'GA_n':          len(ga_costs),
                        'SSM_SA_RPD':    rpd(ssm_sa_t, bks),
                        'MILP_RPD':      rpd(milp_t, bks),
                    })

    if missing_table:
        print(f"  WARNING: {len(missing_table)} keys not found in tables: {missing_table[:5]} …")

    # ── Write CSV ──────────────────────────────────────────────────────────────
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(records[0].keys())
    with open(OUT_CSV, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(records)
    print(f"Written: {OUT_CSV} ({len(records)} rows)")

    # ── Print summary ──────────────────────────────────────────────────────────
    algo_cols   = ['SA_avg_RPD', 'GRASP_avg_RPD', 'GA_avg_RPD', 'SSM_SA_RPD', 'MILP_RPD']
    algo_labels = ['SA',         'GRASP',          'GA',          'SSM-SA',     'MILP']
    print(f"\n{'Algorithm':10s} {'n':>4s} {'mean':>8s} {'std':>8s} {'median':>8s} {'min':>8s} {'max':>8s}")
    print('-' * 58)
    for col, label in zip(algo_cols, algo_labels):
        s = summarise(records, col)
        if s:
            print(f"{label:10s} {s['n']:>4d} "
                  f"{s['mean']:>8.4f} {s['std']:>8.4f} "
                  f"{s['median']:>8.4f} {s['min']:>8.4f} {s['max']:>8.4f}")

    # ── Write supplementary markdown ──────────────────────────────────────────
    write_markdown(records, algo_cols, algo_labels)
    print(f"Written: {OUT_MD}")


def write_markdown(records, algo_cols, algo_labels):
    lines = []
    lines += [
        "# Supplementary Material: RPD Results per Instance",
        "",
        "This file accompanies the paper "
        "*Scheduling rubber shoe sole production on a parallel machine with synchronized interruptions*.",
        "",
        "**RPD** (Relative Percentage Deviation) is defined as:",
        "",
        "```",
        "RPD = (cost − BKS) / BKS × 100",
        "```",
        "",
        "where **BKS** (Best Known Solution) is the minimum feasible cost achieved by any method "
        "across all runs for that instance.",
        "",
        "- `—` for SSM-SA: instance not executed (computationally infeasible for repeated runs).",
        "- `—` for MILP: Gurobi did not find a proven optimal solution within the time limit "
        "(lower bound reported in the paper; feasible value excluded from analysis).",
        "",
    ]

    # Summary table
    lines += [
        "## Summary Statistics (120 instances, both β values)",
        "",
        "| Algorithm | n | Mean RPD | Std RPD | Median RPD | Min RPD | Max RPD |",
        "|-----------|--:|--------:|--------:|-----------:|--------:|--------:|",
    ]
    for col, label in zip(algo_cols, algo_labels):
        s = summarise(records, col)
        if s:
            lines.append(
                f"| {label} | {s['n']} "
                f"| {fmt(s['mean'])} | {fmt(s['std'])} "
                f"| {fmt(s['median'])} | {fmt(s['min'])} | {fmt(s['max'])} |"
            )
    lines.append("")

    # Per-instance table
    lines += [
        "## Per-Instance Average RPD",
        "",
        "| Instance | β | BKS | SA avg | SA std | GRASP avg | GRASP std "
        "| GA avg | GA std | SSM-SA | MILP |",
        "|----------|:-:|----:|-------:|-------:|----------:|----------:"
        "|-------:|-------:|-------:|-----:|",
    ]
    for r in records:
        lines.append(
            f"| {r['instance']} | {r['beta']} | {r['BKS']:.0f} "
            f"| {fmt(r['SA_avg_RPD'])} | {fmt(r['SA_std_RPD'])} "
            f"| {fmt(r['GRASP_avg_RPD'])} | {fmt(r['GRASP_std_RPD'])} "
            f"| {fmt(r['GA_avg_RPD'])} | {fmt(r['GA_std_RPD'])} "
            f"| {fmt(r['SSM_SA_RPD'])} | {fmt(r['MILP_RPD'])} |"
        )
    lines.append("")

    OUT_MD.write_text('\n'.join(lines), encoding='utf-8')


if __name__ == '__main__':
    main()
