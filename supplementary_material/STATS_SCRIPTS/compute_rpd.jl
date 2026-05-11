#!/usr/bin/env julia
"""
compute_rpd.jl — Compute per-run RPD for SA, GRASP, GA and single-value RPD for SSM-SA.

Julia port of compute_rpd.py.

MILP and SSM-SA costs are read from the LaTeX tables in main.tex:
  - Plain number          → optimal / single-run cost
  - "B xxx (F yyy)"      → MILP did not find optimal; excluded from BKS (marked missing)
  - "**"                 → SSM-SA not run for this instance (marked missing)

BKS per instance = min over all available feasible costs.

Required packages (install once):
    import Pkg; Pkg.add(["CSV", "DataFrames"])

Usage:
    julia compute_rpd.jl

Outputs:
    <script_dir>/rpd_results.csv      — per-instance RPD stats (used by statistical_analysis.jl)
    <repo_root>/supplementary_rpd.md  — GitHub supplementary material
"""

using CSV
using DataFrames
using Statistics
using Printf

# ── Paths ──────────────────────────────────────────────────────────────────────
const SCRIPT_DIR = dirname(abspath(@__FILE__))
const REPO_ROOT  = abspath(joinpath(SCRIPT_DIR, "..", ".."))
const AUX_DIR    = joinpath(REPO_ROOT, "aux_temp")
const MAIN_TEX   = joinpath(REPO_ROOT, "main.tex")
const OUT_CSV    = joinpath(SCRIPT_DIR, "rpd_results.csv")
const OUT_MD     = joinpath(REPO_ROOT, "supplementary_rpd.md")

# ── Regex patterns ─────────────────────────────────────────────────────────────
const RE_TABLE_LABEL = r"\\label\{table:combined_methods(_beta6)?_(\d)\}"
const RE_TABLE_BLOCK = r"\\begin\{table\}.*?\\end\{table\}"s
const RE_ORDER       = r"\\multirow\{4\}\{\*\}\{\$O_(\d)\$\}"
const RE_BF          = r"B\s*[\d.]+\s*\(F\s*[\d.]+\)"
const RE_STRUCT      = r"\\(?:toprule|bottomrule|midrule|cmidrule\s*\([^)]*\)\s*\{[^}]*\})"

# ── Cell parsing ───────────────────────────────────────────────────────────────
"""
Parse a cost cell from the LaTeX table.
Returns Float64, or `missing` for non-optimal MILP (B/F format), missing entries (**), or empty.
"""
function parse_t_cell(raw::AbstractString)::Union{Float64,Missing}
    c = strip(raw)
    (isempty(c) || c in ("**", "*", "--")) && return missing
    occursin(RE_BF, c) && return missing
    v = tryparse(Float64, c)
    v === nothing ? missing : v
end

# ── Parse LaTeX tables ─────────────────────────────────────────────────────────
"""
Parse main.tex and extract MILP and SSM-SA costs from result tables.
Returns Dict mapping (scenario, beta, order, p) → NamedTuple(milp, ssm_sa).
"""
function parse_latex_tables(tex_path::String)
    tex    = read(tex_path, String)
    result = Dict{NTuple{4,Int}, @NamedTuple{milp::Union{Float64,Missing}, ssm_sa::Union{Float64,Missing}}}()

    for block_m in eachmatch(RE_TABLE_BLOCK, tex)
        block = block_m.match

        # Only process result-summary tables
        label_m = match(RE_TABLE_LABEL, block)
        label_m === nothing && continue

        beta     = label_m[1] !== nothing ? 6 : 3
        scenario = parse(Int, label_m[2])

        # Extract tabular content
        tab_m = match(r"\\begin\{tabular\}.*?\\end\{tabular\}"s, block)
        tab_m === nothing && continue

        tabular = tab_m.match
        rows    = split(tabular, "\\\\")

        current_order = nothing
        for row in rows
            # Detect order number from \multirow{4}{*}{$O_n$}
            order_m = match(RE_ORDER, row)
            if order_m !== nothing
                current_order = parse(Int, order_m[1])
            end

            # Strip LaTeX structural commands, then split on &
            cleaned = replace(row, RE_STRUCT => "")
            cells   = split(cleaned, "&")
            length(cells) < 4 && continue

            p = tryparse(Int, strip(cells[2]))
            (p === nothing || p ∉ (2, 3, 4, 5)) && continue
            current_order === nothing && continue

            milp_t   = parse_t_cell(cells[3])
            ssm_sa_t = parse_t_cell(cells[4])

            result[(scenario, beta, current_order, p)] = (milp=milp_t, ssm_sa=ssm_sa_t)
        end
    end

    return result
end

# ── Load heuristic run CSVs ────────────────────────────────────────────────────
"""
Load objective-function costs from a heuristic run CSV file.
Returns Vector{Float64} (empty if file is missing or has no 'cost' column).
"""
function load_csv_costs(algo::String, beta::Int, scenario::Int, order::Int, p::Int)::Vector{Float64}
    folder_map = Dict(
        "sa"    => joinpath(AUX_DIR, "logs_sa_beta$(beta)",    "sa_runs"),
        "grasp" => joinpath(AUX_DIR, "logs_grasp_beta$(beta)", "grasp_runs"),
        "ga"    => joinpath(AUX_DIR, "logs_ga_beta$(beta)",    "ga_runs"),
    )
    fpath = joinpath(folder_map[algo], "H_O$(order)_#$(scenario)_$(p)p.csv")
    if !isfile(fpath)
        @warn "Missing CSV: $fpath"
        return Float64[]
    end
    df = CSV.read(fpath, DataFrame; silencewarnings=true)
    "cost" ∉ names(df) && return Float64[]
    return Float64[x for x in df.cost if !ismissing(x)]
end

# ── Statistics helpers ─────────────────────────────────────────────────────────
function rpd(cost::Union{Float64,Missing}, bks::Float64)::Union{Float64,Missing}
    ismissing(cost) && return missing
    bks == 0.0 && return missing
    return (cost - bks) / bks * 100.0
end

function safe_mean(v::Vector{Float64})::Union{Float64,Missing}
    isempty(v) ? missing : mean(v)
end

function safe_std(v::Vector{Float64})::Union{Float64,Missing}
    length(v) < 2 ? missing : std(v)
end

function safe_median(v::Vector{Float64})::Union{Float64,Missing}
    isempty(v) ? missing : median(v)
end

fmt(::Missing, digits::Int=4) = "—"
fmt(v::Real,   digits::Int=4) = @sprintf("%.*f", digits, Float64(v))

function summarise(vals::AbstractVector)
    v = Float64[x for x in vals if !ismissing(x)]
    isempty(v) && return nothing
    return (n=length(v), mean=mean(v), std=std(v), median=median(v), min=minimum(v), max=maximum(v))
end

# ── Write supplementary markdown ───────────────────────────────────────────────
function write_markdown(records::Vector{<:NamedTuple}, algo_cols, algo_labels)
    lines = String[
        "# Supplementary Material: RPD Results per Instance",
        "",
        "This file accompanies the paper " *
        "*Scheduling rubber shoe sole production on a parallel machine with synchronized interruptions*.",
        "",
        "**RPD** (Relative Percentage Deviation) is defined as:",
        "",
        "```",
        "RPD = (cost − BKS) / BKS × 100",
        "```",
        "",
        "where **BKS** (Best Known Solution) is the minimum feasible cost achieved by any method " *
        "across all runs for that instance.",
        "",
        "- `—` for SSM-SA: instance not executed (computationally infeasible for repeated runs).",
        "- `—` for MILP: Gurobi did not find a proven optimal solution within the time limit " *
        "(lower bound reported in the paper; feasible value excluded from analysis).",
        "",
    ]

    # Summary table
    push!(lines, "## Summary Statistics (120 instances, both β values)", "")
    push!(lines,
        "| Algorithm | n | Mean RPD | Std RPD | Median RPD | Min RPD | Max RPD |",
        "|-----------|--:|--------:|--------:|-----------:|--------:|--------:|",
    )
    for (col, label) in zip(algo_cols, algo_labels)
        vals = getfield.(records, col)
        s = summarise(vals)
        s !== nothing && push!(lines,
            "| $label | $(s.n) | $(fmt(s.mean)) | $(fmt(s.std)) | $(fmt(s.median)) | $(fmt(s.min)) | $(fmt(s.max)) |"
        )
    end
    push!(lines, "")

    # Per-instance table
    push!(lines, "## Per-Instance Average RPD", "")
    push!(lines,
        "| Instance | β | BKS | SA avg | SA std | GRASP avg | GRASP std | GA avg | GA std | SSM-SA | MILP |",
        "|----------|:-:|----:|-------:|-------:|----------:|----------:|-------:|-------:|-------:|-----:|",
    )
    for r in records
        push!(lines,
            "| $(r.instance) | $(r.beta) | $(round(Int, r.BKS)) " *
            "| $(fmt(r.SA_avg_RPD)) | $(fmt(r.SA_std_RPD)) " *
            "| $(fmt(r.GRASP_avg_RPD)) | $(fmt(r.GRASP_std_RPD)) " *
            "| $(fmt(r.GA_avg_RPD)) | $(fmt(r.GA_std_RPD)) " *
            "| $(fmt(r.SSM_SA_RPD)) | $(fmt(r.MILP_RPD)) |"
        )
    end
    push!(lines, "")

    write(OUT_MD, join(lines, "\n"))
end

# ── Main ───────────────────────────────────────────────────────────────────────
function main()
    println("Parsing LaTeX tables …")
    table_values = parse_latex_tables(MAIN_TEX)
    println("  Found $(length(table_values)) instance entries in tables.")

    RType = @NamedTuple{
        instance::String, beta::Int, BKS::Float64,
        SA_avg_RPD::Union{Float64,Missing}, SA_std_RPD::Union{Float64,Missing}, SA_n::Int,
        GRASP_avg_RPD::Union{Float64,Missing}, GRASP_std_RPD::Union{Float64,Missing}, GRASP_n::Int,
        GA_avg_RPD::Union{Float64,Missing}, GA_std_RPD::Union{Float64,Missing}, GA_n::Int,
        SSM_SA_RPD::Union{Float64,Missing}, MILP_RPD::Union{Float64,Missing}
    }
    records      = RType[]
    missing_keys = NTuple{4,Int}[]

    for beta in (3, 6), scenario in 1:5, order in 1:3, p in (2, 3, 4, 5)
        key = (scenario, beta, order, p)
        tv  = get(table_values, key, nothing)
        if tv === nothing
            push!(missing_keys, key)
            tv = (milp=missing, ssm_sa=missing)
        end

        sa_costs    = load_csv_costs("sa",    beta, scenario, order, p)
        grasp_costs = load_csv_costs("grasp", beta, scenario, order, p)
        ga_costs    = load_csv_costs("ga",    beta, scenario, order, p)

        all_costs = Float64[]
        ismissing(tv.milp)   || push!(all_costs, tv.milp)
        ismissing(tv.ssm_sa) || push!(all_costs, tv.ssm_sa)
        append!(all_costs, sa_costs)
        append!(all_costs, grasp_costs)
        append!(all_costs, ga_costs)

        if isempty(all_costs)
            @warn "No cost data for $key, skipping."
            continue
        end

        bks = minimum(all_costs)

        sa_rpds    = Float64[rpd(c, bks)::Float64 for c in sa_costs]
        grasp_rpds = Float64[rpd(c, bks)::Float64 for c in grasp_costs]
        ga_rpds    = Float64[rpd(c, bks)::Float64 for c in ga_costs]

        push!(records, (
            instance      = "H_O$(order)_#$(scenario)_$(p)p",
            beta          = beta,
            BKS           = bks,
            SA_avg_RPD    = safe_mean(sa_rpds),
            SA_std_RPD    = safe_std(sa_rpds),
            SA_n          = length(sa_costs),
            GRASP_avg_RPD = safe_mean(grasp_rpds),
            GRASP_std_RPD = safe_std(grasp_rpds),
            GRASP_n       = length(grasp_costs),
            GA_avg_RPD    = safe_mean(ga_rpds),
            GA_std_RPD    = safe_std(ga_rpds),
            GA_n          = length(ga_costs),
            SSM_SA_RPD    = rpd(tv.ssm_sa, bks),
            MILP_RPD      = rpd(tv.milp, bks),
        ))
    end

    if !isempty(missing_keys)
        @warn "$(length(missing_keys)) keys not found in tables: $(missing_keys[1:min(5,end)]) …"
    end

    # ── Write CSV ────────────────────────────────────────────────────────────
    mkpath(dirname(OUT_CSV))
    df = DataFrame(records)
    CSV.write(OUT_CSV, df; missingstring="—")
    println("Written: $OUT_CSV ($(nrow(df)) rows)")

    # ── Print summary ─────────────────────────────────────────────────────────
    algo_cols   = (:SA_avg_RPD, :GRASP_avg_RPD, :GA_avg_RPD, :SSM_SA_RPD, :MILP_RPD)
    algo_labels = ("SA", "GRASP", "GA", "SSM-SA", "MILP")

    @printf("\n%-10s %4s %8s %8s %8s %8s %8s\n", "Algorithm", "n", "mean", "std", "median", "min", "max")
    println(repeat('-', 58))
    for (col, label) in zip(algo_cols, algo_labels)
        s = summarise(getfield.(records, col))
        s !== nothing && @printf("%-10s %4d %8.4f %8.4f %8.4f %8.4f %8.4f\n",
            label, s.n, s.mean, s.std, s.median, s.min, s.max)
    end

    # ── Write supplementary markdown ──────────────────────────────────────────
    write_markdown(records, algo_cols, algo_labels)
    println("Written: $OUT_MD")
end

main()
