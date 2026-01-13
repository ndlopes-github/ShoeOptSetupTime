#= Copyright (C) 2024
Nuno David Lopes.
Created:  2025/11/21
Last changed - N. Lopes: 2026/01/13 11:53:31

Run script for Greedy algorithm
This script executes a deterministic greedy constructive heuristic for quick baseline solutions.
The greedy algorithm provides fast, repeatable results for single-mold scenarios.

USAGE:
    julia --project scripts/run_greedy.jl [OPTIONS]

OPTIONS:
    --file=PATH          Load custom instance file from data/settings/
    --help, -h           Display this help message

EXAMPLES:
    julia --project scripts/run_greedy.jl                      # Run default instance
    julia --project scripts/run_greedy.jl --file=E_O1_#1_2p.jl  # Run custom instance

DESCRIPTION:
    The greedy algorithm is a fast deterministic constructive heuristic that:
    - Always selects the largest remaining subjob
    - Places it on the shelf with minimum resulting load
    - Provides a quick baseline solution for comparison with other methods
    
RESTRICTIONS:
    Only supports single-mold scenarios (all jobs must have o[j] == 1).
    For multi-mold instances, use GRASP or other algorithms instead.
=#

using DrWatson
@quickactivate "ShoeOptSetupTime"
using Printf

"""
    parse_greedy_arguments()::Union{String, Nothing}

Parse command-line arguments for the Greedy run script.

# Supported Arguments
- `--file=PATH`: Load custom instance file from data/settings/ directory
- `--help, -h`: Display help message and exit

# Returns
- `Union{String, Nothing}`: Custom file path if specified, else nothing

# Default Behavior
If no arguments provided, loads default instance: H_O2_#2_3p.jl
"""
function parse_greedy_arguments()
    custom_file = nothing

    for arg in ARGS
        if startswith(arg, "--file=")
            custom_file = String(split(arg, "=")[2])
        elseif arg == "--help" || arg == "-h"
            println("Greedy Algorithm Run Script")
            println("=" ^ 80)
            println("USAGE:")
            println("    julia --project scripts/run_greedy.jl [OPTIONS]")
            println()
            println("OPTIONS:")
            println("    --file=PATH          Load custom instance file from data/settings/")
            println("    --help, -h           Display this help message")
            println()
            println("EXAMPLES:")
            println("    julia --project scripts/run_greedy.jl")
            println("    julia --project scripts/run_greedy.jl --file=E_O1_#1_2p.jl")
            println()
            println("DESCRIPTION:")
            println("    Fast deterministic greedy heuristic for baseline solutions.")
            println("    Only supports single-mold scenarios (all o[j] == 1).")
            println("=" ^ 80)
            exit(0)
        else
            @warn "Unknown argument: $(arg). Use --help for usage information."
        end
    end

    return custom_file
end

# Parse arguments
custom_file = parse_greedy_arguments()

# Load settings FIRST, before including the module
# The instance file must define `order_dict` with required keys:
# - `:g` - Job identifiers
# - `:n` - Job quantities  
# - `:o` - Molds per job (must all be 1 for greedy)
# - `:p` - Number of parallel shelves/machines
# - `:α` - Shelf penalty coefficient
# - `:β` - Job delay penalty coefficient
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

# Validate single-mold requirement
if any(x -> x != 1, order_dict[:o])
    @error "Greedy algorithm only supports single-mold scenarios (all o[j] == 1)"
    @error "Current instance has multi-mold jobs. Use GRASP or other algorithms instead."
    exit(1)
end

# Now include and use the GRASP module's greedy wrapper
include(scriptsdir("grasp.jl"))
using .Grasp

println("=" ^ 80)
println("Greedy Algorithm (Deterministic)")
println("=" ^ 80)

# Run Greedy and measure time
elapsed_time = @elapsed begin
    result = Grasp.run_greedy(order_dict)
end

println("\n" * "=" ^ 80)
println("GREEDY EXECUTION COMPLETED")
println("=" ^ 80)
println(@sprintf("Cost: %.6f", result.cost))
println(@sprintf("Shelf levels (m): %d", result.m))
println(@sprintf("Maximum shelf load: %.6f", result.max_sum))
println(@sprintf("Total execution time: %.3fs", elapsed_time))
println("=" ^ 80)
