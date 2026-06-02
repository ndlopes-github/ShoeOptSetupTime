#!/usr/bin/env julia
"""
compute_rpd.jl — Compute per-run RPD for SA, GRASP, GA and single-value RPD for SSM-SA.

MILP and SSM-SA costs are read from:
  STATS_SCRIPTS/tables_04_to_13_combined.csv
  (columns: beta, scenario, order, p, t_MILP, t_SSM_SA, ...)
  "**" in t_MILP or t_SSM_SA → instance not solved / excluded from BKS.

Per-run SA / GRASP / GA costs are read from CSV files produced by the batch
runs and stored inside the supplementary_material archives:
  supplementary_material/SA/sa_beta{3,6}_logs.tar.bz2
    └── logs_sa_beta{N}/sa_runs/H_O{order}_#{scenario}_{p}p.csv
  supplementary_material/GRASP/grasp_beta{3,6}_logs.tar.bz2
    └── logs_grasp_beta{N}/grasp_runs/H_O{order}_#{scenario}_{p}p.csv
  supplementary_material/GA/ga_beta{3,6}_logs.tar.bz2
    └── logs_ga_beta{N}/ga_runs/H_O{order}_#{scenario}_{p}p.csv

Extract the archives to a directory of your choice and set RUN_DIR below to
that directory, so that the following paths exist:
  <RUN_DIR>/logs_sa_beta3/sa_runs/H_O1_#1_2p.csv  (etc.)

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
const REPO_ROOT = abspath(joinpath(SCRIPT_DIR, "..", ".."))
const TABLES_CSV = joinpath(SCRIPT_DIR, "tables_04_to_13_combined.csv")
const OUT_CSV = joinpath(SCRIPT_DIR, "rpd_results.csv")
const OUT_MD = joinpath(REPO_ROOT, "supplementary_rpd.md")

# ── Directory where the heuristic run archives were extracted ──────────────────
# Extract the six archives in supplementary_material/{SA,GRASP,GA}/ to a
# directory and set RUN_DIR to that path.  Example:
#   tar -xjf supplementary_material/SA/sa_beta3_logs.tar.bz2 -C /tmp/runs
#   tar -xjf supplementary_material/SA/sa_beta6_logs.tar.bz2 -C /tmp/runs
#   ... (repeat for GRASP and GA archives)
#   Then set: RUN_DIR = "/tmp/runs"
const RUN_DIR = joinpath(REPO_ROOT, "supplementary_material")   # adjust if needed

# ── Load MILP / SSM-SA values from combined CSV ────────────────────────────────
"""
Parse tables_04_to_13_combined.csv.
Returns Dict mapping (scenario, beta, order, p) → NamedTuple(milp, ssm_sa).
"**" cells → missing.
"""
function load_table_values(csv_path::String)
    df = CSV.read(csv_path, DataFrame; missingstring=["", "**", "*", "--"])
    result = Dict{NTuple{4,Int},@NamedTuple{milp::Union{Float64,Missing}, ssm_sa::Union{Float64,Missing}}}()
    for row in eachrow(df)
        key = (Int(row.scenario), Int(row.beta), Int(row.order), Int(row.p))
        milp = ismissing(row.t_MILP) ? missing : tryparse(Float64, string(row.t_MILP))
        ssm_sa = ismissing(row.t_SSM_SA) ? missing : tryparse(Float64, string(row.t_SSM_SA))
        result[key] = (milp=something(milp, missing), ssm_sa=something(ssm_sa, missing))
    end
    return result
end

# ── Load heuristic run CSVs ────────────────────────────────────────────────────
"""
Load objective-function costs from a heuristic run CSV file.
Returns Vector{Float64} (empty if file is missing or has no 'cost' column).
If the expected directory does not exist at all, prints a one-time hint about
the archives that need to be extracted.
"""
const _MISSING_HINT_SHOWN = Ref(false)

function load_csv_costs(algo::String, beta::Int, scenario::Int, order::Int, p::Int)::Vector{Float64}
    run_subdir = Dict(
        "sa" => joinpath(RUN_DIR, "logs_sa_beta$(beta)", "sa_runs"),
        "grasp" => joinpath(RUN_DIR, "logs_grasp_beta$(beta)", "grasp_runs"),
        "ga" => joinpath(RUN_DIR, "logs_ga_beta$(beta)", "ga_runs"),
    )
    dir = run_subdir[algo]
    fpath = joinpath(dir, "H_O$(order)_#$(scenario)_$(p)p.csv")

    if !isdir(dir) && !_MISSING_HINT_SHOWN[]
        _MISSING_HINT_SHOWN[] = true
        println("""
\n  ╔══════════════════════════════════════════════════════════════════╗
  ║  Per-run CSV files not found.                                    ║
  ║                                                                  ║
  ║  The run data is compressed in:                                  ║
  ║    supplementary_material/SA/sa_beta3_logs.tar.bz2              ║
  ║    supplementary_material/SA/sa_beta6_logs.tar.bz2              ║
  ║    supplementary_material/GRASP/grasp_beta3_logs.tar.bz2        ║
  ║    supplementary_material/GRASP/grasp_beta6_logs.tar.bz2        ║
  ║    supplementary_material/GA/ga_beta3_logs.tar.bz2              ║
  ║    supplementary_material/GA/ga_beta6_logs.tar.bz2              ║
  ║                                                                  ║
  ║  Extract them all to a directory, e.g.:                          ║
  ║    for f in supplementary_material/{SA,GRASP,GA}/*.tar.bz2; do  ║
  ║        tar -xjf "\$f" -C /tmp/runs                               ║
  ║    done                                                          ║
  ║                                                                  ║
  ║  Then set RUN_DIR at the top of this script to that path,        ║
  ║  e.g.:  const RUN_DIR = "/tmp/runs"                              ║
  ╚══════════════════════════════════════════════════════════════════╝
""")
    end

    if !isfile(fpath)
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
fmt(v::Real, digits::Int=4) = @sprintf("%.*f", digits, Float64(v))

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
        "This file accompanies the paper "*"*Scheduling rubber shoe sole production on a parallel machine with synchronized interruptions*.",
        "",
        "**RPD** (Relative Percentage Deviation) is defined as:",
        "",
        "```",
        "RPD = (cost − BKS) / BKS × 100",
        "```",
        "",
        "where **BKS** (Best Known Solution) is the minimum feasible cost achieved by any method "*"across all runs for that instance.",
        "",
        "- `—` for SSM-SA: heuristic not applicable — no feasible two-subset job partition exists for the instance (denoted `**` in the main paper).",
        "- `—` for MILP: Gurobi did not find a proven optimal solution within the time limit "*"(lower bound reported in the paper; feasible value excluded from analysis).",
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
    println("Loading table values from $(basename(TABLES_CSV)) …")
    table_values = load_table_values(TABLES_CSV)
    println("  Found $(length(table_values)) instance entries.")

    RType = @NamedTuple{
        instance::String, beta::Int, BKS::Float64,
        SA_avg_RPD::Union{Float64,Missing}, SA_std_RPD::Union{Float64,Missing}, SA_n::Int,
        GRASP_avg_RPD::Union{Float64,Missing}, GRASP_std_RPD::Union{Float64,Missing}, GRASP_n::Int,
        GA_avg_RPD::Union{Float64,Missing}, GA_std_RPD::Union{Float64,Missing}, GA_n::Int,
        SSM_SA_RPD::Union{Float64,Missing}, MILP_RPD::Union{Float64,Missing}
    }
    records = RType[]
    missing_keys = NTuple{4,Int}[]

    for beta in (3, 6), scenario in 1:5, order in 1:3, p in (2, 3, 4, 5)
        key = (scenario, beta, order, p)
        tv = get(table_values, key, nothing)
        if tv === nothing
            push!(missing_keys, key)
            tv = (milp=missing, ssm_sa=missing)
        end

        sa_costs = load_csv_costs("sa", beta, scenario, order, p)
        grasp_costs = load_csv_costs("grasp", beta, scenario, order, p)
        ga_costs = load_csv_costs("ga", beta, scenario, order, p)

        all_costs = Float64[]
        ismissing(tv.milp) || push!(all_costs, tv.milp)
        ismissing(tv.ssm_sa) || push!(all_costs, tv.ssm_sa)
        append!(all_costs, sa_costs)
        append!(all_costs, grasp_costs)
        append!(all_costs, ga_costs)

        if isempty(all_costs)
            @warn "No cost data for $key, skipping."
            continue
        end

        bks = minimum(all_costs)

        sa_rpds = Float64[rpd(c, bks)::Float64 for c in sa_costs]
        grasp_rpds = Float64[rpd(c, bks)::Float64 for c in grasp_costs]
        ga_rpds = Float64[rpd(c, bks)::Float64 for c in ga_costs]

        push!(records, (
            instance="H_O$(order)_#$(scenario)_$(p)p",
            beta=beta,
            BKS=bks,
            SA_avg_RPD=safe_mean(sa_rpds),
            SA_std_RPD=safe_std(sa_rpds),
            SA_n=length(sa_costs),
            GRASP_avg_RPD=safe_mean(grasp_rpds),
            GRASP_std_RPD=safe_std(grasp_rpds),
            GRASP_n=length(grasp_costs),
            GA_avg_RPD=safe_mean(ga_rpds),
            GA_std_RPD=safe_std(ga_rpds),
            GA_n=length(ga_costs),
            SSM_SA_RPD=rpd(tv.ssm_sa, bks),
            MILP_RPD=rpd(tv.milp, bks),
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
    algo_cols = (:SA_avg_RPD, :GRASP_avg_RPD, :GA_avg_RPD, :SSM_SA_RPD, :MILP_RPD)
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
