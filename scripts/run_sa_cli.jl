#= Copyright (C) 2024
Nuno David Lopes.
Created:  2026/04/20
Last changed - N. Lopes:2026/04/28 19:46:59

CLI runner for pure Simulated Annealing.
Loads a settings instance, overrides SA parameters with irace best-configuration
tuned values, and runs SimulatedAnnealing independently Nruns times, reporting
the best solution found. Per-run costs are saved to a CSV file under
data/exp_pro/sa_runs/ for statistical analysis.

Usage:
    julia --project scripts/run_sa_cli.jl --instance=H_O1_#1_2p.jl
    julia --project scripts/run_sa_cli.jl --instance H_O1_#1_2p.jl
=#
using DrWatson
@quickactivate "SoftIdea"

using Printf, Random, Logging, CSV, DataFrames

# ─── Irace best-configuration tuned parameters ───────────────────────────────
# Beta = 3 Irace Parameters
#const IRACE_T0 = 48.2242
#const IRACE_Tf = 0.0791
#const IRACE_Nit = 1995
#const IRACE_Tj = 13
#const IRACE_NRUNS = 181


# Beta = 6 Irace parameters
const IRACE_T0 = 0.7724
const IRACE_Tf = 0.1615
const IRACE_Nit = 1885
const IRACE_Tj = 11
const IRACE_NRUNS = 196

# ─── Output configuration ────────────────────────────────────────────────────
const OUTPUT_DIR = datadir("exp_pro", "sa_runs")

"""
Parse CLI arguments and return the instance settings filename.

# Arguments
- `args`: Command-line arguments vector

# Returns
- `String`: Instance settings filename
"""
function parse_instance_arg(args::Vector{String})::String
    positional_args = String[]

    for (index, arg) in pairs(args)
        if arg == "--instance"
            index == lastindex(args) && error("Missing value for --instance: EXIT")
            return strip(args[index+1])
        end

        if startswith(arg, "--instance=")
            value = strip(split(arg, "=", limit=2)[2])
            isempty(value) && error("Missing value for --instance: EXIT")
            return value
        end

        startswith(arg, "--") || push!(positional_args, strip(arg))
    end

    isempty(positional_args) && error(
        "No instance argument provided. Use --instance <file>, --instance=<file>, or a positional file name: EXIT",
    )

    return first(positional_args)
end

"""Parse optional --beta argument. Returns the value as a Float64, or nothing."""
function parse_beta_arg(args::Vector{String})
    for (index, arg) in pairs(args)
        if arg == "--beta"
            index == lastindex(args) && error("Missing value for --beta: EXIT")
            return parse(Float64, args[index+1])
        end
        if startswith(arg, "--beta=")
            return parse(Float64, split(arg, "=", limit=2)[2])
        end
    end
    return nothing
end

instance = parse_instance_arg(ARGS)
beta_override = parse_beta_arg(ARGS)
settings_path = datadir("settings", instance)
isfile(settings_path) || error("Settings file not found: $(settings_path): EXIT")

@info "Running pure SA for instance: $instance"
@info "Parameters (irace best config): T0=$(IRACE_T0)  Tf=$(IRACE_Tf)  Nit=$(IRACE_Nit)  Tj=$(IRACE_Tj)  Nruns=$(IRACE_NRUNS)"
!isnothing(beta_override) && @info "Beta override: β = $beta_override (instance value bypassed)"

include(settings_path)
# order_dict is now defined in the current scope via include

# ─── Override with irace-tuned parameters ────────────────────────────────────
order_dict[:T0] = IRACE_T0
order_dict[:Tf] = IRACE_Tf
order_dict[:Nit] = IRACE_Nit
order_dict[:Tj] = IRACE_Tj
!isnothing(beta_override) && (order_dict[:β] = beta_override)

include(scriptsdir("simulated_annealing.jl"))
using .SimulatedAnnealing

@unpack α, β = order_dict

# ─── Warm-up run (discard compilation overhead) ──────────────────────────────
with_logger(NullLogger()) do
    SimulatedAnnealing.run(order_dict; log_every=0)
end

# ─── Independent runs ────────────────────────────────────────────────────────
best_overall_shelves = nothing
best_overall_cost = Inf
best_overall_m = 0
best_overall_max_sum = 0.0
best_run_idx = 0
total_time = 0.0
run_times = Float64[]
run_costs = Float64[]
run_ms = Int[]

println("="^80)
println("Instance  : $instance")
println("Parameters: T0=$(IRACE_T0)  Tf=$(IRACE_Tf)  Nit=$(IRACE_Nit)  Tj=$(IRACE_Tj)  Nruns=$(IRACE_NRUNS)")
println("="^80)

for run_idx in 1:IRACE_NRUNS
    run_time = @elapsed begin
        shelves = with_logger(NullLogger()) do
            SimulatedAnnealing.run(order_dict; log_every=0)
        end
        max_sum, m, cost = SimulatedAnnealing.get_cost_from_shelves(shelves, α, β)

        push!(run_costs, cost)
        push!(run_ms, m)
        if cost < best_overall_cost
            global best_overall_cost = cost
            global best_overall_m = m
            global best_overall_max_sum = max_sum
            global best_overall_shelves = deepcopy(shelves)
            global best_run_idx = run_idx
        end
    end

    push!(run_times, run_time)
    global total_time += run_time

    if run_idx % 10 == 0 || run_idx == IRACE_NRUNS
        @sprintf("%3d/%d  best_so_far=%.6f  elapsed=%.1fs",
            run_idx, IRACE_NRUNS, best_overall_cost, total_time) |> println
    end
end

# ─── Final report ─────────────────────────────────────────────────────────────
println("\n" * "="^80)
println("RESULTS: $instance")
println("="^80)
println("Best cost      : $best_overall_cost")
println("Best m         : $best_overall_m")
println("Best max_sum   : $best_overall_max_sum")
println("Best run index : $best_run_idx / $(IRACE_NRUNS)")
println(@sprintf("Total time     : %.2fs  |  avg %.2fs/run  |  min %.2fs  |  max %.2fs",
    total_time,
    total_time / IRACE_NRUNS,
    minimum(run_times),
    maximum(run_times)))
println("="^80)

# ─── Save per-run costs ──────────────────────────────────────────────────────
mkpath(OUTPUT_DIR)
csv_path = joinpath(OUTPUT_DIR, replace(instance, ".jl" => ".csv"))
CSV.write(csv_path, DataFrame(run_idx=1:length(run_costs), cost=run_costs, m=run_ms, time=run_times))
println(@sprintf("Per-run costs  : saved to %s  (%d rows)", csv_path, length(run_costs)))
