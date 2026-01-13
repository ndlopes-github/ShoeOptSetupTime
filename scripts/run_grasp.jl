#= Copyright (C) 2024
Nuno David Lopes.
Created:  2025/11/25

Run script for GRASP algorithm
This script executes the Greedy Randomized Adaptive Search Procedure (GRASP) on a specified instance.
The instance is configured via a settings file.

USAGE:
    julia --project scripts/run_grasp.jl [OPTIONS]

OPTIONS:
    --file=PATH          Load custom instance file from data/settings/
    --iterations=N       Override number of GRASP iterations (default: from settings file)
    --help, -h           Display this help message

EXAMPLES:
    julia --project scripts/run_grasp.jl                           # Run default instance
    julia --project scripts/run_grasp.jl --file=E_O1_#1_2p.jl      # Run custom instance
    julia --project scripts/run_grasp.jl --iterations=50           # Override iterations

DESCRIPTION:
    GRASP is a multi-start metaheuristic that combines greedy construction with randomization.
    Each iteration builds a solution using randomized greedy choices, and the best solution
    across all iterations is returned. This implementation is tailored to the parallel machine
    scheduling problem with simultaneous interruptions.
=#

using DrWatson
@quickactivate "ShoeOptSetupTime"
using Printf

"""
    parse_grasp_arguments()::Tuple

Parse command-line arguments for the GRASP run script.

# Supported Arguments
- `--file=PATH`: Load custom instance file from data/settings/ directory
- `--iterations=N`: Override the number of GRASP iterations from settings file
- `--help, -h`: Display help message and exit

# Returns
- `(custom_file, iterations_override)::Tuple`:
  - `custom_file::Union{String, Nothing}`: Custom file path if specified, else nothing
  - `iterations_override::Union{Int, Nothing}`: Iteration count override or nothing

# Default Behavior
If no arguments provided, loads default instance: H_O2_#2_3p.jl with Nit from settings file
"""
function parse_grasp_arguments()
    custom_file = nothing
    iterations_override = nothing

    for arg in ARGS
        if startswith(arg, "--file=")
            custom_file = String(split(arg, "=")[2])
        elseif startswith(arg, "--iterations=")
            iterations_override = parse(Int, split(arg, "=")[2])
        elseif arg == "--help" || arg == "-h"
            println("GRASP Algorithm Run Script")
            println("=" ^ 80)
            println("USAGE:")
            println("    julia --project scripts/run_grasp.jl [OPTIONS]")
            println()
            println("OPTIONS:")
            println("    --file=PATH          Load custom instance file from data/settings/")
            println("    --iterations=N       Override number of GRASP iterations")
            println("    --help, -h           Display this help message")
            println()
            println("EXAMPLES:")
            println("    julia --project scripts/run_grasp.jl")
            println("    julia --project scripts/run_grasp.jl --file=E_O1_#1_2p.jl")
            println("    julia --project scripts/run_grasp.jl --iterations=50")
            println("=" ^ 80)
            exit(0)
        else
            @warn "Unknown argument: $(arg). Use --help for usage information."
        end
    end

    return custom_file, iterations_override
end

# Parse arguments
custom_file, iterations_override = parse_grasp_arguments()

# Load settings
# The instance file must define `order_dict` with required keys:
# - `:g` - Job identifiers
# - `:n` - Job quantities
# - `:o` - Molds per job
# - `:p` - Number of parallel shelves/machines
# - `:α` - Shelf penalty coefficient
# - `:β` - Job delay penalty coefficient
# - `:Nit` - Number of GRASP iterations
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

# Now include and use the GRASP module
include(scriptsdir("grasp.jl"))
using .Grasp

println("=" ^ 80)
println("GRASP Algorithm")
println("=" ^ 80)
if iterations_override !== nothing
    println("Iterations: $(iterations_override) (override from command line)")
else
    println("Iterations: $(order_dict[:Nit]) (from settings file)")
end
println("=" ^ 80)

# Run GRASP with timing
elapsed_time = @elapsed begin
    if iterations_override !== nothing
        result = Grasp.run(order_dict; iterations=iterations_override)
    else
        result = Grasp.run(order_dict)
    end
end

println("\n" * "=" ^ 80)
println("GRASP EXECUTION COMPLETED")
println("=" ^ 80)
println("Total execution time: $(round(elapsed_time, digits=3))s")
println("=" ^ 80)
