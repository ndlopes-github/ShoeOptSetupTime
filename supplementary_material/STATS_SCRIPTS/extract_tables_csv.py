#!/usr/bin/env python3
"""
extract_tables_csv.py — Extract raw t / s / time data from tables 4–13 in main.tex.

Tables 4–8  : beta=3, scenarios 1–5 (label: table:combined_methods_{1..5})
Tables 9–13 : beta=6, scenarios 1–5 (label: table:combined_methods_beta6_{1..5})

Outputs (in the same folder as this script):
  table04_scenario1_beta3.csv  …  table13_scenario5_beta6.csv   (one per table)
  tables_04_to_13_combined.csv                                   (all 120 rows)

Column layout per table:
  Scenario 1 (has Greedy G): 6 methods × 3 metric groups  (t / s / time)
  Scenarios 2–5:             5 methods × 3 metric groups
"""

import csv
import re
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT  = SCRIPT_DIR.parent.parent
MAIN_TEX   = REPO_ROOT / "main.tex"

# ── Regex ──────────────────────────────────────────────────────────────────────
RE_LABEL   = re.compile(r'\\label\{table:combined_methods(_beta6)?_(\d)\}')
RE_ORDER   = re.compile(r'\$O_(\d)\$')
RE_STRUCT  = re.compile(r'\\(?:toprule|bottomrule|midrule|hline|cmidrule\s*\([^)]*\)\s*\{[^}]*\})')
RE_STRIP   = re.compile(r'\\(?:multirow|multicolumn)\{[^}]*\}\{[^}]*\}\{[^}]*\}')
RE_TABLE   = re.compile(r'\\begin\{table\}.*?\\end\{table\}', re.DOTALL)

# ── Cell cleaning ──────────────────────────────────────────────────────────────
def clean(cell: str) -> str:
    """Strip LaTeX markup from a cell value; return plain text."""
    c = cell.strip()
    # Remove \new{...} wrapper
    c = re.sub(r'\\new\{([^}]*)\}', r'\1', c)
    # Remove color/text commands
    c = re.sub(r'\\(?:color|textbf|textit|emph)\{[^}]*\}\{([^}]*)\}', r'\1', c)
    # Remove remaining braces from single-arg commands
    c = re.sub(r'\\[a-zA-Z]+\{([^}]*)\}', r'\1', c)
    c = re.sub(r'[{}]', '', c)
    return c.strip()


def parse_value(raw: str) -> str:
    """
    Return a clean string representation:
      plain number  → the number as string
      B x (F y)     → "B x (F y)"   (MILP non-optimal)
      ** or *       → same
      empty         → ""
    """
    c = clean(raw)
    if not c:
        return ""
    # Normalise B x (F y) format
    if re.search(r'B\s*[\d.]+\s*\(F', c):
        nums = re.findall(r'[\d.]+', c)
        return f"B {nums[0]} (F {nums[1]})" if len(nums) >= 2 else c
    if c in ('**', '*'):
        return c
    # Try numeric
    try:
        v = float(c)
        return str(int(v)) if v == int(v) else str(v)
    except ValueError:
        return c


# ── Main parser ────────────────────────────────────────────────────────────────
def parse_tables(tex_path: Path):
    """
    Returns list of dicts, one per data row across all 10 tables.
    Dict keys: table_num, beta, scenario, order, p,
               t_MILP, t_SSM_SA, t_SA, t_GRASP, t_GA, [t_G],
               s_MILP, s_SSM_SA, s_SA, s_GRASP, s_GA, [s_G],
               time_MILP, time_SSM_SA, time_SA, time_GRASP, time_GA, [time_G]
    """
    tex     = tex_path.read_text(encoding='utf-8')
    records = []
    table_num = 3  # will increment to 4 at first match

    for block_m in re.finditer(RE_TABLE, tex):
        block     = block_m.group(0)
        label_m   = RE_LABEL.search(block)
        if not label_m:
            continue

        table_num += 1
        beta       = 6 if label_m.group(1) else 3
        scenario   = int(label_m.group(2))
        has_G      = (scenario == 1)
        n_methods  = 6 if has_G else 5

        # Extract tabular body
        tab_m = re.search(r'\\begin\{tabular\}.*?\\end\{tabular\}', block, re.DOTALL)
        if not tab_m:
            continue

        tabular = tab_m.group(0)
        rows    = tabular.split('\\\\')

        current_order = None

        for row in rows:
            # Detect order from $O_n$
            order_m = RE_ORDER.search(row)
            if order_m:
                current_order = int(order_m.group(1))

            # Strip structural commands
            cleaned = RE_STRUCT.sub('', row)
            cells   = cleaned.split('&')

            # Need: 1 (order/empty) + 1 (p) + n_methods×3 (t+s+time)
            expected = 2 + n_methods * 3
            if len(cells) < expected:
                continue

            p_val = parse_value(cells[1])
            try:
                p = int(p_val)
            except ValueError:
                continue
            if p not in (2, 3, 4, 5) or current_order is None:
                continue

            # Slice value cells (after order and p)
            vals = [parse_value(cells[i]) for i in range(2, 2 + n_methods * 3)]

            # t block: vals[0..n_methods-1]
            # s block: vals[n_methods..2*n_methods-1]
            # time block: vals[2*n_methods..3*n_methods-1]
            t    = vals[:n_methods]
            s    = vals[n_methods:2 * n_methods]
            time = vals[2 * n_methods:3 * n_methods]

            methods = ['MILP', 'SSM_SA', 'SA', 'GRASP', 'GA'] + (['G'] if has_G else [])

            rec = {
                'table_num': table_num,
                'beta':      beta,
                'scenario':  scenario,
                'order':     current_order,
                'p':         p,
            }
            for idx, m in enumerate(methods):
                rec[f't_{m}']    = t[idx]
                rec[f's_{m}']    = s[idx]
                rec[f'time_{m}'] = time[idx]

            records.append(rec)

    return records


# ── Write CSVs ─────────────────────────────────────────────────────────────────
def write_csvs(records):
    # --- individual per-table CSVs ---
    from itertools import groupby

    key_fn = lambda r: (r['table_num'], r['beta'], r['scenario'])
    by_table = {}
    for r in records:
        k = key_fn(r)
        by_table.setdefault(k, []).append(r)

    table_names = {}  # table_num -> filename
    for (tnum, beta, scenario), rows in sorted(by_table.items()):
        fname = SCRIPT_DIR / f'table{tnum:02d}_scenario{scenario}_beta{beta}.csv'
        table_names[tnum] = fname.name
        fieldnames = list(rows[0].keys())
        with open(fname, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(rows)
        print(f"Written: {fname} ({len(rows)} rows)")

    # --- combined CSV ---
    if not records:
        return
    combined = SCRIPT_DIR / 'tables_04_to_13_combined.csv'
    # Use the superset of all keys (scenario 1 has G columns, others don't)
    all_keys = []
    seen = set()
    for r in records:
        for k in r:
            if k not in seen:
                all_keys.append(k)
                seen.add(k)
    with open(combined, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=all_keys, extrasaction='ignore', restval='')
        w.writeheader()
        w.writerows(records)
    print(f"\nWritten: {combined} ({len(records)} total rows)")


# ── Report timing stats ────────────────────────────────────────────────────────
def report_timing_stats(records):
    """Print min/max/mean for each heuristic's time column."""
    import math

    def stats(vals):
        v = [float(x) for x in vals if x and x not in ('*', '**') and
             not x.startswith('B')]
        if not v:
            return None
        return dict(n=len(v), min=min(v), max=max(v),
                    mean=sum(v) / len(v))

    print('\n=== Timing statistics (seconds) ===')
    methods = ['SA', 'GRASP', 'GA', 'G']
    header  = f"{'Method':10s} {'n':>4s} {'min':>8s} {'max':>8s} {'mean':>8s}"
    print(header)
    print('-' * len(header))
    for m in methods:
        key  = f'time_{m}'
        vals = [r.get(key, '') for r in records]
        s    = stats(vals)
        if s:
            print(f"{m:10s} {s['n']:>4d} {s['min']:>8.3f} {s['max']:>8.3f} {s['mean']:>8.3f}")


def main():
    print(f"Parsing {MAIN_TEX} …")
    records = parse_tables(MAIN_TEX)
    print(f"Found {len(records)} data rows across tables 4–13.\n")

    write_csvs(records)
    report_timing_stats(records)


if __name__ == '__main__':
    main()
