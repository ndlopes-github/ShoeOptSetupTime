#= Copyright (C) 2025
Nuno David Lopes.
Created:  2026/04/28

CLI runner for Genetic Algorithm.
Loads a settings instance, overrides GA parameters with best-configuration
values, and runs GeneticAlgorithm independently Nruns times, reporting the
best solution found.

Usage:
    julia --project scripts/run_ga_cli.jl --instance=H_O1_#1_2p.jl
    julia --project scripts/run_ga_cli.jl --instance H_O1_#1_2p.jl
=#
using DrWatson
@quickactivate "SoftIdea"

using Printf, Random, Logging, CSV, DataFrames

# ─── Output configuration ─────────────────────────────────────────────
const OUTPUT_DIR = datadir("exp_pro", "ga_runs")

# ─── Best-configuration parameters ──────────────────────────────────────────
# NOTE: GA has not yet been through irace tuning. These are working defaults.
const GA_POP_SIZE = 50
const GA_CLONE_THRESHOLD = 0.1
const GA_Nit = 100   # overrides order_dict[:Nit] (generations)
const GA_NRUNS = 100

"""
Parse CLI arguments and return the instance settings filename.
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

@info "Running GA for instance: $instance"
@info "Parameters: pop_size=$(GA_POP_SIZE)  clone_threshold=$(GA_CLONE_THRESHOLD)  Nit=$(GA_Nit)  Nruns=$(GA_NRUNS)"
!isnothing(beta_override) && @info "Beta override: β = $beta_override (instance value bypassed)"

include(settings_path)
# order_dict is now defined in the current scope via include

# ─── Override with configured parameters ─────────────────────────────────────
order_dict[:Nit] = GA_Nit
!isnothing(beta_override) && (order_dict[:β] = beta_override)

include(scriptsdir("simulated_annealing.jl"))
using .SimulatedAnnealing: get_cost_from_shelves

include(scriptsdir("genetic_algorithm.jl"))
using .GeneticAlgorithm

@unpack α, β = order_dict

# ─── Warm-up run (discard JIT compilation overhead) ──────────────────────────
with_logger(NullLogger()) do
    GeneticAlgorithm.run_ga(order_dict;
        pop_size=GA_POP_SIZE,
        clone_threshold=GA_CLONE_THRESHOLD, log_every=0)
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
println("Parameters: pop_size=$(GA_POP_SIZE)  clone_threshold=$(GA_CLONE_THRESHOLD)  Nit=$(GA_Nit)  Nruns=$(GA_NRUNS)")
println("="^80)

for run_idx in 1:GA_NRUNS
    run_time = @elapsed begin
        shelves = with_logger(NullLogger()) do
            GeneticAlgorithm.run_ga(order_dict;
                pop_size=GA_POP_SIZE,
                clone_threshold=GA_CLONE_THRESHOLD, log_every=0)
        end
        max_sum, m, cost = get_cost_from_shelves(shelves, α, β)

        push!(run_costs, Float64(cost))
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

    if run_idx % 5 == 0 || run_idx == GA_NRUNS
        @sprintf("%3d/%d  best_so_far=%.6f  elapsed=%.1fs",
            run_idx, GA_NRUNS, best_overall_cost, total_time) |> println
    end
end

# ─── Final report ─────────────────────────────────────────────────────────────
println("\n" * "="^80)
println("RESULTS: $instance")
println("="^80)
println("Best cost      : $best_overall_cost")
println("Best m         : $best_overall_m")
println("Best max_sum   : $best_overall_max_sum")
println("Best run index : $best_run_idx / $(GA_NRUNS)")
println(@sprintf("Total time     : %.2fs  |  avg %.2fs/run  |  min %.2fs  |  max %.2fs",
    total_time,
    total_time / GA_NRUNS,
    minimum(run_times),
    maximum(run_times)))
println("="^80)

# ─── Save per-run costs ───────────────────────────────────────────────────────
mkpath(OUTPUT_DIR)
csv_path = joinpath(OUTPUT_DIR, replace(instance, ".jl" => ".csv"))
CSV.write(csv_path, DataFrame(run_idx=1:length(run_costs), cost=run_costs, m=run_ms, time=run_times))
println(@sprintf("Per-run costs  : saved to %s  (%d rows)", csv_path, length(run_costs)))
