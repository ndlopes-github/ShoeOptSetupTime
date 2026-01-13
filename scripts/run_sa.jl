#= Copyright (C) 2024
Nuno David Lopes.
Created:  2024/04/09
Last changed - N. Lopes: 2026/01/13 12:15:05

Run script for Simulated Annealing algorithm
This script executes multiple independent runs of Simulated Annealing (SA) and reports
the best solution found. Multiple runs help assess solution quality and algorithm consistency.

USAGE:
    julia --project scripts/run_sa.jl [OPTIONS]

OPTIONS:
    --file=PATH          Load custom instance file from data/settings/
    --runs=N             Number of independent SA runs (default: 100)
    --log-every=N        Logging frequency within each SA run (default: 10, 0 to disable)
    --help, -h           Display this help message

EXAMPLES:
    julia --project scripts/run_sa.jl                         # Run default (100 runs)
    julia --project scripts/run_sa.jl --runs=50               # Run with 50 trials
    julia --project scripts/run_sa.jl --file=E_O1_#1_2p.jl    # Custom instance
    julia --project scripts/run_sa.jl --runs=10 --log-every=0 # 10 runs, silent mode

DESCRIPTION:
    Simulated Annealing is a probabilistic metaheuristic that explores the solution space
    by accepting both improving and worsening moves with temperature-dependent probability.
    Multiple independent runs are used to:
    - Find better solutions through stochastic exploration
    - Assess solution quality variability
    - Provide statistical evidence of performance
=#

using DrWatson
@quickactivate "ShoeOptSetupTime"
using Printf

"""
    parse_sa_arguments()::Tuple

Parse command-line arguments for the Simulated Annealing run script.

# Supported Arguments
- `--file=PATH`: Load custom instance file from data/settings/ directory
- `--runs=N`: Number of independent SA runs (default: 100)
- `--log-every=N`: Logging frequency within each SA run (default: 10)
- `--help, -h`: Display help message and exit

# Returns
- `(custom_file, num_runs, log_every)::Tuple`:
  - `custom_file::Union{String, Nothing}`: Custom file path if specified, else nothing
  - `num_runs::Int`: Number of independent runs (default: 100)
  - `log_every::Int`: Logging frequency (default: 10)

# Default Behavior
If no arguments provided:
- Loads default instance: H_O2_#2_3p.jl
- Runs 100 independent trials
- Logs progress every 10 iterations within each run
"""
function parse_sa_arguments()
    custom_file = nothing
    num_runs = 100  # default
    log_every = 10  # default

    for arg in ARGS
        if startswith(arg, "--file=")
            custom_file = String(split(arg, "=")[2])
        elseif startswith(arg, "--runs=")
            num_runs = parse(Int, split(arg, "=")[2])
            if num_runs <= 0
                @error "Number of runs must be positive"
                exit(1)
            end
        elseif startswith(arg, "--log-every=")
            log_every = parse(Int, split(arg, "=")[2])
        elseif arg == "--help" || arg == "-h"
            println("Simulated Annealing Run Script")
            println("=" ^ 80)
            println("USAGE:")
            println("    julia --project scripts/run_sa.jl [OPTIONS]")
            println()
            println("OPTIONS:")
            println("    --file=PATH          Load custom instance file from data/settings/")
            println("    --runs=N             Number of independent SA runs (default: 100)")
            println("    --log-every=N        Logging frequency per run (default: 10, 0=silent)")
            println("    --help, -h           Display this help message")
            println()
            println("EXAMPLES:")
            println("    julia --project scripts/run_sa.jl")
            println("    julia --project scripts/run_sa.jl --runs=50")
            println("    julia --project scripts/run_sa.jl --file=E_O1_#1_2p.jl --runs=10")
            println()
            println("DESCRIPTION:")
            println("    Runs multiple independent Simulated Annealing trials and reports")
            println("    the best solution found across all runs.")
            println("=" ^ 80)
            exit(0)
        else
            @warn "Unknown argument: $(arg). Use --help for usage information."
        end
    end

    return custom_file, num_runs, log_every
end

# Parse arguments
custom_file, num_runs, log_every = parse_sa_arguments()

# Load settings
# The instance file must define `order_dict` with required keys:
# - `:g` - Job identifiers
# - `:n` - Job quantities
# - `:o` - Molds per job
# - `:p` - Number of parallel shelves/machines
# - `:α` - Shelf penalty coefficient
# - `:β` - Job delay penalty coefficient
# - `:T0` - Initial temperature (SA parameter)
# - `:αT` - Temperature cooling rate (SA parameter)
# - `:maxIt` - Maximum iterations per SA run
if custom_file !== nothing
    # Load custom instance file
    instance_file = datadir("settings", custom_file)
    if !isfile(instance_file)
        @error "Custom instance file not found: $(instance_file)"
        exit(1)
    end
    println("Loading custom instance: $(custom_file)")
    include(instance_file)
else
    # Load default instance
    include(datadir("settings", "H_O2_#2_3p.jl"))
    println("Loading instance: H_O2_#2_3p.jl (default)")
end

# Now include and use the Simulated Annealing module
include(scriptsdir("simulated_annealing.jl"))
using .SimulatedAnnealing

@unpack α, β = order_dict

println("=" ^ 80)
println("Simulated Annealing")
println("=" ^ 80)
println("Number of independent runs: $(num_runs)")
println("Logging frequency: $(log_every > 0 ? "every $(log_every) iterations" : "disabled")")
println("=" ^ 80)

best_overall_shelves = nothing
best_overall_cost = Inf
best_overall_m = 0
best_overall_max_sum = 0
best_run_idx = 0
total_time = 0.0
run_times = Float64[]

for run_idx in 1:num_runs
    println("\n" * "=" ^ 80)
    println("Run $(run_idx)/$(num_runs)")
    println("=" ^ 80)
    
    # Run SA with timing
    run_time = @elapsed begin
        shelves = SimulatedAnnealing.run(order_dict; log_every=log_every)
    end
    push!(run_times, run_time)
    global total_time += run_time
    
    # Compute cost for this run
    max_sum, m, cost = SimulatedAnnealing.get_cost_from_shelves(shelves, α, β)
    
    println("\nRun $(run_idx) results: cost=$(cost), m=$(m), max_sum=$(max_sum), time=$(round(run_time, digits=3))s")
    
    # Update best solution if this run is better
    if cost < best_overall_cost
        global best_overall_cost = cost
        global best_overall_m = m
        global best_overall_max_sum = max_sum
        global best_overall_shelves = deepcopy(shelves)
        global best_run_idx = run_idx
        println("  ✓ New best solution found!")
    end
end

# Report final results
println("\n" * "=" ^ 80)
println("FINAL RESULTS AFTER $(num_runs) RUNS")
println("=" ^ 80)
println("Best solution found in run $(best_run_idx)")
println("Best cost: $(best_overall_cost)")
println("Number of non-empty slots (m): $(best_overall_m)")
println("Maximum cumulative quantity: $(best_overall_max_sum)")
println("\nTiming Statistics:")
println("  Total time: $(round(total_time, digits=3))s")
println("  Average time per run: $(round(total_time/num_runs, digits=3))s")
println("  Min time: $(round(minimum(run_times), digits=3))s")
println("  Max time: $(round(maximum(run_times), digits=3))s")
println("\nBest partition:")
for (i, shelf) in enumerate(best_overall_shelves)
    if !isempty(shelf)
        println("  Shelf $i: ", join(["(job=$(j.jobid), mold=$(j.moldid), qty=$(j.quantity))" for j in shelf], ", "))
    end
end
println("=" ^ 80)
