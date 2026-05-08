#=  Copyright (C) 2025
Nuno David Lopes.
Created: 2025/12/11
=#

using DrWatson
@quickactivate "SoftIdea"

using Logging, LoggingExtras, Printf

# Include the GA module
include(scriptsdir("genetic_algorithm.jl"))
using .GeneticAlgorithm

include(scriptsdir("simulated_annealing.jl"))
using .SimulatedAnnealing: get_cost_from_shelves

# --- Configuration ---
# Beta override: set to nothing to use the value from the instance file,
# or set a number (e.g. 12) to bypass it.
const BETA_OVERRIDE = nothing

# Log file setup
log_path = projectdir("logs", "run_ga.log")
mkpath(dirname(log_path))
const TEE_LOGGER = TeeLogger(
    MinLevelLogger(FileLogger(log_path), Logging.Info),
    MinLevelLogger(ConsoleLogger(stderr), Logging.Debug)
)
global_logger(TEE_LOGGER)

# --- Select and Load Instance ---
# You can change this to any instance file from `data/settings`
instance_file = "H_O3_#3_3p.jl"
instance_path = datadir("settings", instance_file)

if !isfile(instance_path)
    @error "Instance file not found: $instance_path"
else
    @info "Loading instance: $instance_file"
    include(instance_path)
end

if !isnothing(BETA_OVERRIDE)
    order_dict[:β] = BETA_OVERRIDE
    @info "Beta override active: β = $BETA_OVERRIDE (instance value bypassed)"
end

# --- Main Execution ---
function main(order_dict)
    # order_dict is passed as a parameter from the top level
    # Nit from order_dict is used as the number of GA generations.

    # --- GA Parameters ---
    # These can be tuned for better performance
    pop_size = 400
    clone_threshold = 0.1
    log_every = 20

    # --- Warm-up run (discard JIT compilation overhead) ---
    with_logger(NullLogger()) do
        run_ga(order_dict; pop_size=pop_size,
            clone_threshold=clone_threshold, log_every=0)
    end

    # --- Run the Algorithm ---
    elapsed_time = @elapsed begin
        best_shelves = run_ga(
            order_dict;
            pop_size=pop_size,
            clone_threshold=clone_threshold,
            log_every=log_every
        )
    end

    @unpack α, β = order_dict
    best_max_sum, best_m, best_cost = get_cost_from_shelves(best_shelves, α, β)

    # --- Display Results ---
    println("\n" * "="^80)
    println("RESULTS: $instance_file")
    println("="^80)
    println("Best cost      : $best_cost")
    println("Best m         : $best_m")
    println("Best max_sum   : $best_max_sum")
    println(@sprintf("Elapsed time   : %.2fs", elapsed_time))
    println("Best shelf assignment:")
    for (i, shelf) in enumerate(best_shelves)
        if !isempty(shelf)
            println("  Shelf $i: ", join(["(job=$(j.jobid), mold=$(j.moldid), qty=$(j.quantity))" for j in shelf], ", "))
        end
    end
    println("="^80)

    # You can also save the results if needed, for example:
    # results = @dict instance_file best_cost best_shelves elapsed_time
    # safesave(datadir("sims", "ga_results.jld2"), results)
end

# --- Run ---
if abspath(PROGRAM_FILE) == @__FILE__
    main(order_dict)
end
