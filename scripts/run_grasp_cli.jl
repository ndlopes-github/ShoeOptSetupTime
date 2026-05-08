#= Copyright (C) 2024
Nuno David Lopes.
Created:  2026/04/27

CLI runner for GRASP.
Loads a settings instance, overrides Nit with the irace best-configuration
tuned value, and runs Grasp.run once (which performs Nit internal constructions
and returns the best solution found). Per-iteration costs are saved to a CSV
file under data/exp_pro/grasp_runs/ for statistical analysis.

Usage:
    julia --project scripts/run_grasp_cli.jl --instance=H_O1_#1_2p.jl
    julia --project scripts/run_grasp_cli.jl --instance H_O1_#1_2p.jl
=#
using DrWatson
@quickactivate "SoftIdea"

using Printf, Logging, CSV, DataFrames

# ─── Irace best-configuration tuned parameters ───────────────────────────────
# Beta = 3 Irace Parameters
const IRACE_Nit = 144

# ─── Output configuration ────────────────────────────────────────────────────
const OUTPUT_DIR = datadir("exp_pro", "grasp_runs")

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

@info "Running GRASP for instance: $instance"
@info "Parameters (irace best config): Nit=$(IRACE_Nit)"
!isnothing(beta_override) && @info "Beta override: β = $beta_override (instance value bypassed)"

include(settings_path)
# order_dict is now defined in the current scope via include

# ─── Override with irace-tuned parameters ────────────────────────────────────
order_dict[:Nit] = IRACE_Nit
!isnothing(beta_override) && (order_dict[:β] = beta_override)

include(scriptsdir("simulated_annealing.jl"))
using .SimulatedAnnealing: get_cost_from_shelves

include(scriptsdir("grasp.jl"))
using .Grasp

@unpack α, β = order_dict

# ─── Warm-up run (discard compilation overhead) ──────────────────────────────
with_logger(NullLogger()) do
    Grasp.run(order_dict; log_every=0)
end

# ─── Main run ────────────────────────────────────────────────────────────────
println("="^80)
println("Instance  : $instance")
println("Parameters: Nit=$(IRACE_Nit)")
println("="^80)

elapsed_time = @elapsed begin
    grasp_result = with_logger(NullLogger()) do
        Grasp.run(order_dict; log_every=0, return_costs=true)
    end
end

best_shelves = grasp_result.shelves
run_costs = grasp_result.run_costs
run_ms = grasp_result.run_ms
run_times = grasp_result.run_times
best_max_sum, best_m, best_cost = get_cost_from_shelves(best_shelves, α, β)

# ─── Final report ─────────────────────────────────────────────────────────────
println("\n" * "="^80)
println("RESULTS: $instance")
println("="^80)
println("Best cost      : $best_cost")
println("Best m         : $best_m")
println("Best max_sum   : $best_max_sum")
println(@sprintf("Elapsed time   : %.2fs", elapsed_time))
println("="^80)

# ─── Save per-iteration costs ────────────────────────────────────────────────
mkpath(OUTPUT_DIR)
csv_path = joinpath(OUTPUT_DIR, replace(instance, ".jl" => ".csv"))
CSV.write(csv_path, DataFrame(iter_idx=1:length(run_costs), cost=run_costs, m=run_ms, time=run_times))
println(@sprintf("Per-iter costs : saved to %s  (%d rows)", csv_path, length(run_costs)))
