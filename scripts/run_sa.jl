#= Copyright (C) 2024
Nuno David Lopes.
Created:  2024/04/09
Last changed - N. Lopes: 2026/01/12 15:52:06
=#

using DrWatson
@quickactivate "ShoeOptSetupTime"
using Printf

# Load settings
include(datadir("settings", "H_O2_#2_3p.jl"))

# Now include and use the (renamed) heuristic module
include(scriptsdir("simulated_annealing.jl"))
using .SimulatedAnnealing

# Configuration: number of independent runs
const NUM_RUNS = 100

@unpack α, β = order_dict

println("=" ^ 80)
println("Running Simulated Annealing with $(NUM_RUNS) independent runs")
println("=" ^ 80)

best_overall_shelves = nothing
best_overall_cost = Inf
best_overall_m = 0
best_overall_max_sum = 0
best_run_idx = 0
total_time = 0.0
run_times = Float64[]


for run_idx in 1:NUM_RUNS
    println("\n" * "=" ^ 80)
    println("Run $(run_idx)/$(NUM_RUNS)")
    println("=" ^ 80)
    
    # Run SA with timing
    run_time = @elapsed begin
        shelves = SimulatedAnnealing.run(order_dict; log_every=10)
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
println("FINAL RESULTS AFTER $(NUM_RUNS) RUNS")
println("=" ^ 80)
println("Best solution found in run $(best_run_idx)")
println("Best cost: $(best_overall_cost)")
println("Number of non-empty slots (m): $(best_overall_m)")
println("Maximum cumulative quantity: $(best_overall_max_sum)")
println("\nTiming Statistics:")
println("  Total time: $(round(total_time, digits=3))s")
println("  Average time per run: $(round(total_time/NUM_RUNS, digits=3))s")
println("  Min time: $(round(minimum(run_times), digits=3))s")
println("  Max time: $(round(maximum(run_times), digits=3))s")
println("\nBest partition:")
for (i, shelf) in enumerate(best_overall_shelves)
    if !isempty(shelf)
        println("  Shelf $i: ", join(["(job=$(j.jobid), mold=$(j.moldid), qty=$(j.quantity))" for j in shelf], ", "))
    end
end
println("=" ^ 80)
