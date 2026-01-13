#= Copyright (C) 2024
Nuno David Lopes.
Created:  2025/11/21
Last changed - N. Lopes: 2026/01/13 11:42:53

Run script for Split-Solve-Merge MILP algorithm
This script executes the Split-Solve-Merge heuristic on a specified instance.
The instance is configured via a settings file (E_* for exact, H_* for heuristic instances).

USAGE:
    julia --project scripts/run_ssm_milp.jl [OPTIONS]

OPTIONS:
    exact           Load exact problem instance (E_O2_#2_3p.jl)
    heuristic       Load heuristic problem instance (H_O2_#2_3p.jl) [default]
    --help, -h      Display this help message
    --file=PATH     Load custom instance file from data/settings/

EXAMPLES:
    julia --project scripts/run_ssm_milp.jl              # Run default heuristic instance
    julia --project scripts/run_ssm_milp.jl exact        # Run exact instance
    julia --project scripts/run_ssm_milp.jl --file=E_O1_#1_2p.jl

DESCRIPTION:
    The Split-Solve-Merge algorithm is a novel heuristic designed for handling
    large-scale instances of the parallel machine scheduling problem with
    simultaneous interruptions. It partitions the problem, solves sub-problems
    independently, and merges the solutions.
=#

using DrWatson
@quickactivate "ShoeOptSetupTime"
using Printf

"""
    parse_ssm_arguments()::Tuple

Parse command-line arguments for the Split-Solve-Merge MILP run script.

# Supported Arguments
- `exact`: Load exact problem instance (E_O2_#2_3p.jl)
- `heuristic`: Load heuristic instance (H_O2_#2_3p.jl) [default]
- `--file=PATH`: Load custom instance file from data/settings/ directory
- `--help, -h`: Display help message and exit

# Returns
- `(instance_type, custom_file)::Tuple`:
  - `instance_type::String`: Either "exact" or "heuristic"
  - `custom_file::Union{String, Nothing}`: Custom file path if specified, else nothing

# Default Behavior
If no arguments provided, loads default heuristic instance: H_O2_#2_3p.jl
"""
function parse_ssm_arguments()
    instance_type = "heuristic"  # default
    custom_file = nothing

    for arg in ARGS
        if arg == "exact"
            instance_type = "exact"
        elseif arg == "heuristic"
            instance_type = "heuristic"
        elseif startswith(arg, "--file=")
            custom_file = String(split(arg, "=")[2])
        elseif arg == "--help" || arg == "-h"
            println("Split-Solve-Merge MILP Run Script")
            println("=" ^ 80)
            println("USAGE:")
            println("    julia --project scripts/run_ssm_milp.jl [OPTIONS]")
            println()
            println("OPTIONS:")
            println("    exact               Load exact problem instance (E_O2_#2_3p.jl)")
            println("    heuristic           Load heuristic problem instance (H_O2_#2_3p.jl) [default]")
            println("    --file=PATH         Load custom instance file from data/settings/")
            println("    --help, -h          Display this help message")
            println()
            println("EXAMPLES:")
            println("    julia --project scripts/run_ssm_milp.jl")
            println("    julia --project scripts/run_ssm_milp.jl exact")
            println("    julia --project scripts/run_ssm_milp.jl --file=E_O1_#1_2p.jl")
            println("=" ^ 80)
            exit(0)
        else
            @warn "Unknown argument: $(arg). Use --help for usage information."
        end
    end

    return instance_type, custom_file
end

# Parse arguments
instance_type, custom_file = parse_ssm_arguments()

# Load settings FIRST, before including the module
# The instance file must define `order_dict` with required keys:
# - `:Oid`: Order identifier
# - `:order`: List of job orders
# - `:α`: Shelf weight penalty
# - `:β`: Job delay penalty

if custom_file !== nothing
    # Load custom instance file
    instance_file = datadir("settings", custom_file)
    if !isfile(instance_file)
        @error "Custom instance file not found: $(instance_file)"
        exit(1)
    end
    println("Loading custom instance: $(custom_file)")
    include(instance_file)
elseif instance_type == "exact"
    # Load exact instance
    include(datadir("settings", "E_O2_#2_3p.jl"))
    println("Loading instance: E_O2_#2_3p.jl (exact)")
else
    # Load heuristic instance (default)
    include(datadir("settings", "H_O2_#2_3p.jl"))
    println("Loading instance: H_O2_#2_3p.jl (heuristic)")
end

# Now include and use the module
include(scriptsdir("split_solve_merge_milp.jl"))
using .SplitSolveMergeMILP

println("=" ^ 80)
println("Split-Solve-Merge MILP Algorithm")
println("=" ^ 80)

# Run SSM with timing
elapsed_time = @elapsed begin
    SplitSolveMergeMILP.run(order_dict)
end

println("\n" * "=" ^ 80)
println("SPLIT-SOLVE-MERGE EXECUTION COMPLETED")
println("=" ^ 80)
println("Total execution time: $(round(elapsed_time, digits=3))s")
println("=" ^ 80)
