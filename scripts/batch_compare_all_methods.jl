#= Copyright (C) 2024
Nuno David Lopes.
Created:  2025/11/28
Comprehensive batch script to compare all methods:
1. MILP (E_* files, exact solver)
2. Split-Solve-Merge (H_* files, heuristic with Pg>1)
3. Simulated Annealing (E_* files, 100 runs)
4. GRASP (E_* files, single run)

Output: One CSV table per scenario (#1, #2, etc.)
Output: One CSV table with all scenarios
=#

using DrWatson
@quickactivate "ShoeOptSetupTime"
using Printf
using DataFrames
using CSV
using Statistics
using Dates

# Load algorithm modules
include(scriptsdir("simulated_annealing.jl"))
using .SimulatedAnnealing

include(scriptsdir("grasp.jl"))
using .Grasp

include(scriptsdir("split_solve_merge_milp.jl"))
using .SplitSolveMergeMILP

# Parse command-line arguments
"""
    parse_arguments()::Tuple

Parse command-line arguments for batch comparison script.

# Supported Arguments
- `--limit=N, -l=N`: Process only first N instances (for testing/dry runs)
- `--beta=VALUE, -b=VALUE`: Override beta penalty value from settings files
- `--only-file=PATH`: Run only instances listed in file (Order,Scenario,P format)
- `--skip-milp`: Skip MILP solver (mark as Skipped in output)
- `--skip-ssm`: Skip Split-Solve-Merge heuristic (mark as Skipped in output)
- `--skip-sa`: Skip Simulated Annealing (mark as Skipped in output)
- `--skip-grasp`: Skip GRASP (mark as Skipped in output)
- `--help, -h`: Display help message and exit

# Returns
- `(dry_run_limit, skip_methods, beta_override, only_file)` tuple:
  - `dry_run_limit::Union{Int, Nothing}`: Instance limit or nothing if not set
  - `skip_methods::Set{String}`: Set of method names to skip
  - `beta_override::Union{Float64, Nothing}`: Beta value override or nothing
  - `only_file::Union{String, Nothing}`: Path to whitelist file or nothing

# Examples
```julia
# julia batch_compare_all_methods.jl --limit=5 --skip-milp
# Process only first 5 instances, skipping MILP solver
```
"""
function parse_arguments()
    dry_run_limit = nothing
    skip_methods = Set{String}()
    beta_override = nothing
    only_file = nothing
    
    for arg in ARGS
        if startswith(arg, "--limit=")
            dry_run_limit = parse(Int, split(arg, "=")[2])
        elseif startswith(arg, "-l=")
            dry_run_limit = parse(Int, split(arg, "=")[2])
        elseif startswith(arg, "--beta=")
            beta_override = parse(Float64, split(arg, "=")[2])
        elseif startswith(arg, "-b=")
            beta_override = parse(Float64, split(arg, "=")[2])
        elseif startswith(arg, "--only-file=")
            only_file = String(split(arg, "=")[2])
        elseif arg == "--skip-milp"
            push!(skip_methods, "MILP")
        elseif arg == "--skip-ssm"
            push!(skip_methods, "SSM")
        elseif arg == "--skip-sa"
            push!(skip_methods, "SA")
        elseif arg == "--skip-grasp"
            push!(skip_methods, "GRASP")
        elseif arg == "--help" || arg == "-h"
            println("Usage: julia batch_compare_all_methods.jl [OPTIONS]")
            println("Options:")
            println("  --limit=N, -l=N    Process only first N instances (for testing)")
            println("  --beta=VALUE, -b=VALUE  Override beta valuComplete with a help system if e from settings files")
            println("  --only-file=PATH   Run only instances listed in file (Order,Scenario,P format)")
            println("  --skip-milp        Skip MILP solver (mark as Skipped in output)")
            println("  --skip-ssm         Skip Split-Solve-Merge (mark as Skipped in output)")
            println("  --skip-sa          Skip Simulated Annealing (mark as Skipped in output)")
            println("  --skip-grasp       Skip GRASP (mark as Skipped in output)")
            println("  --help, -h         Show this help message")
            exit(0)
        end
    end
    return dry_run_limit, skip_methods, beta_override, only_file
end

# Configuration
const NUM_SA_RUNS = 100  # Number of independent SA runs per instance
const OUTPUT_DIR = datadir("exp_pro")
const OUTPUT_CSV = joinpath(OUTPUT_DIR, "all_methods_comparison.csv")
const PROGRESS_LOG = joinpath(OUTPUT_DIR, "batch_progress.log")
const ERROR_LOG = joinpath(OUTPUT_DIR, "batch_errors.log")
const DRY_RUN_LIMIT, SKIP_METHODS, BETA_OVERRIDE, ONLY_FILE = parse_arguments()  # Limit instances and methods to skip

"""
    format_time(t::Float64)::String

Format elapsed time for display with adaptive precision based on magnitude.

# Arguments
- `t::Float64`: Time in seconds

# Returns
- `String`: Formatted time with appropriate decimal places

# Examples
```julia
format_time(0.005)   # "0.005"
format_time(0.5)     # "0.50"
format_time(50.5)    # "50.5"
format_time(150.0)   # "150"
```
"""
function format_time(t::Float64)
    if t < 0.01
        return @sprintf("%.3f", t)
    elseif t < 1.0
        return @sprintf("%.2f", t)
    elseif t < 100.0
        return @sprintf("%.1f", t)
    else
        return @sprintf("%.0f", t)
    end
end

"""
    log_error(log_file::String, method::String, instance_id::String, error_obj, stacktrace_info=nothing)

Log detailed error information including timestamp, error type, message, and optional stacktrace.
Appends to existing log file maintaining session history.

# Arguments
- `log_file::String`: Path to error log file
- `method::String`: Name of the algorithm/method that failed
- `instance_id::String`: Unique identifier for the instance being processed
- `error_obj`: The exception object
- `stacktrace_info=nothing`: Optional stacktrace frames for debugging

# Side Effects
- Appends formatted error entry to `log_file`
"""
function log_error(log_file::String, method::String, instance_id::String, error_obj, stacktrace_info=nothing)
    timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    open(log_file, "a") do io
        println(io, "="^80)
        println(io, "[$(timestamp)] ERROR in $(method)")
        println(io, "Instance: $(instance_id)")
        println(io, "Error type: $(typeof(error_obj))")
        println(io, "Error message: $(error_obj)")
        if stacktrace_info !== nothing
            println(io, "Stacktrace:")
            for (i, frame) in enumerate(stacktrace_info)
                println(io, "  $(i). $(frame)")
            end
        end
        println(io, "="^80)
        println(io, "")
    end
end

"""
    parse_instance_name(filename::String)::Union{NamedTuple, Nothing}

Parse instance filename to extract configuration parameters.
Expected format: `[E|H]_O<order>_#<scenario>_<p>p.jl`

# Arguments
- `filename::String`: Instance filename, e.g., "E_O1_#2_2p.jl" or "H_O1_#2_2p.jl"

# Returns
- `NamedTuple` with fields: `(type, order, scenario, p)` if match successful
- `nothing` if filename does not match expected pattern

# Examples
```julia
parse_instance_name("E_O1_#2_2p.jl")
# (type="E", order="O1", scenario=2, p=2)
```
"""
function parse_instance_name(filename::String)
    # E_O1_#2_2p.jl or H_O1_#2_2p.jl -> type=E/H, order=O1, scenario=2, p=2
    m = match(r"([EH])_(O\d+)_#(\d+)_(\d+)p\.jl", filename)
    if m !== nothing
        return (type=m.captures[1], order=m.captures[2], scenario=parse(Int, m.captures[3]), p=parse(Int, m.captures[4]))
    end
    return nothing
end

"""
    run_milp(order_dict::Dict)::NamedTuple

Run the exact MILP solver on a given instance via Split-Solve-Merge module.
Catches and logs errors gracefully.

# Arguments
- `order_dict::Dict`: Instance configuration dictionary with keys including `:Oid`, `:α`, `:β`

# Returns
- `NamedTuple` with fields: `(cost, m, time, error)`
  - `cost`: Objective value (numeric or "Err")
  - `m`: Number of shelf levels (integer or "Err")
  - `time`: Execution time in seconds (numeric or "Err")
  - `error::Bool`: Whether an error occurred

# Side Effects
- Logs errors to ERROR_LOG if solver fails
"""
function run_milp(order_dict::Dict)
    @unpack α, β = order_dict
    
    try
        cost, m, time = SplitSolveMergeMILP.run(order_dict)
        return (cost=cost, m=m, time=time, error=false)
    catch e
        instance_id = "$(order_dict[:Oid])"
        @warn "MILP failed for $(instance_id): $e"
        log_error(ERROR_LOG, "MILP", instance_id, e, catch_backtrace())
        return (cost="Err", m="Err", time="Err", error=true)
    end
end

"""
    run_split_solve_merge(order_dict::Dict)::NamedTuple

Run the Split-Solve-Merge heuristic on a given instance for large-scale problems.
Designed for instances with Pg > 1.

# Arguments
- `order_dict::Dict`: Instance configuration dictionary with keys including `:Oid`, `:α`, `:β`

# Returns
- `NamedTuple` with fields: `(cost, m, time, error)`
  - `cost`: Objective value (numeric or "Err")
  - `m`: Number of shelf levels (integer or "Err")
  - `time`: Execution time in seconds (numeric or "Err")
  - `error::Bool`: Whether an error occurred

# Side Effects
- Logs errors to ERROR_LOG if heuristic fails
"""
function run_split_solve_merge(order_dict::Dict)
    @unpack α, β = order_dict
    
    try
        cost, m, time = SplitSolveMergeMILP.run(order_dict)
        return (cost=cost, m=m, time=time, error=false)
    catch e
        instance_id = "$(order_dict[:Oid])"
        @warn "Split-Solve-Merge failed for $(instance_id): $e"
        log_error(ERROR_LOG, "Split-Solve-Merge", instance_id, e, catch_backtrace())
        return (cost="Err", m="Err", time="Err", error=true)
    end
end

"""
    run_sa_multiple(order_dict::Dict, num_runs::Int)::NamedTuple

Run Simulated Annealing multiple times and track the best solution found.
Useful for assessing consistency and variability of the metaheuristic.

# Arguments
- `order_dict::Dict`: Instance configuration dictionary with keys including `:Oid`, `:α`, `:β`
- `num_runs::Int`: Number of independent SA runs to execute

# Returns
- `NamedTuple` with fields: `(cost, m, time, error)`
  - `cost`: Best objective value across all runs (numeric or "Err")
  - `m`: Shelf levels for best solution (integer or "Err")
  - `time`: Total cumulative time across all runs (numeric or "Err")
  - `error::Bool`: Whether an error occurred

# Side Effects
- Logs errors to ERROR_LOG if any run fails
- Time includes all runs, useful for comparing total computational budget
"""
function run_sa_multiple(order_dict::Dict, num_runs::Int)
    @unpack α, β = order_dict
    
    try
        best_cost = Inf
        best_m = 0
        total_time = 0.0
        
        for run_idx in 1:num_runs
            run_time = @elapsed begin
                shelves = SimulatedAnnealing.run(order_dict; log_every=0)
            end
            total_time += run_time
            
            max_sum, m, cost = SimulatedAnnealing.get_cost_from_shelves(shelves, α, β)
            
            if cost < best_cost
                best_cost = cost
                best_m = m
            end
        end
        
        return (cost=best_cost, m=best_m, time=total_time, error=false)
    catch e
        instance_id = "$(order_dict[:Oid])"
        @warn "SA failed for $(instance_id): $e"
        log_error(ERROR_LOG, "Simulated Annealing", instance_id, e, catch_backtrace())
        return (cost="Err", m="Err", time="Err", error=true)
    end
end

"""
    run_grasp_once(order_dict::Dict)::NamedTuple

Run the Greedy Randomized Adaptive Search Procedure (GRASP) once on a given instance.
Provides a single-run deterministic comparison point.

# Arguments
- `order_dict::Dict`: Instance configuration dictionary with keys including `:Oid`, `:α`, `:β`

# Returns
- `NamedTuple` with fields: `(cost, m, time, error)`
  - `cost`: Objective value (numeric or "Err")
  - `m`: Number of shelf levels (integer or "Err")
  - `time`: Execution time in seconds (numeric or "Err")
  - `error::Bool`: Whether an error occurred

# Side Effects
- Logs errors to ERROR_LOG if GRASP fails
"""
function run_grasp_once(order_dict::Dict)
    @unpack α, β = order_dict
    
    try
        run_time = @elapsed begin
            shelves = Grasp.run(order_dict; log_every=0)
        end
        
        max_sum, m, cost = SimulatedAnnealing.get_cost_from_shelves(shelves, α, β)
        
        return (cost=cost, m=m, time=run_time, error=false)
    catch e
        instance_id = "$(order_dict[:Oid])"
        @warn "GRASP failed for $(instance_id): $e"
        log_error(ERROR_LOG, "GRASP", instance_id, e, catch_backtrace())
        return (cost="Err", m="Err", time="Err", error=true)
    end
end

"""
    load_completed_instances(log_file::String)::Set

Read progress log to identify instances already processed.
Enables resuming batch runs without recomputing completed instances.

# Arguments
- `log_file::String`: Path to progress log file (format: "Order,Scenario,P" per line)

# Returns
- `Set{Tuple{String,Int,Int}}`: Set of (order, scenario, p) tuples representing completed instances
- Empty set if log file does not exist

# Log File Format
Comments start with '#'. Data lines have format: `Order,Scenario,P` (comma-separated)
"""
function load_completed_instances(log_file::String)
    completed = Set{Tuple{String,Int,Int}}()
    if isfile(log_file)
        for line in eachline(log_file)
            line = strip(line)
            if isempty(line) || startswith(line, "#")
                continue
            end
            # Parse format: Order,Scenario,P
            parts = split(line, ",")
            if length(parts) == 3
                try
                    order = String(strip(parts[1]))
                    scenario = parse(Int, strip(parts[2]))
                    p = parse(Int, strip(parts[3]))
                    push!(completed, (order, scenario, p))
                catch e
                    @warn "Could not parse progress log line: $(line)"
                end
            end
        end
        println("  Found $(length(completed)) completed instances in progress log")
    end
    return completed
end

"""
    load_only_instances(only_file::String)::Set

Read a whitelist file to filter which instances should be processed.
Useful for running specific instances or retrying failed cases.

# Arguments
- `only_file::String`: Path to file containing allowed instances (format: "Order,Scenario,P" per line)

# Returns
- `Set{Tuple{String,Int,Int}}`: Set of (order, scenario, p) tuples to process
- Empty set if file does not exist or cannot be read

# Log File Format
Comments start with '#'. Data lines have format: `Order,Scenario,P` (comma-separated)
"""
function load_only_instances(only_file::String)
    only_instances = Set{Tuple{String,Int,Int}}()
    if isfile(only_file)
        for line in eachline(only_file)
            line = strip(line)
            if isempty(line) || startswith(line, "#")
                continue
            end
            # Parse format: Order,Scenario,P
            parts = split(line, ",")
            if length(parts) == 3
                try
                    order = String(strip(parts[1]))
                    scenario = parse(Int, strip(parts[2]))
                    p = parse(Int, strip(parts[3]))
                    push!(only_instances, (order, scenario, p))
                catch e
                    @warn "Could not parse only-file line: $(line)"
                end
            end
        end
        println("  Found $(length(only_instances)) instances to run from only-file")
    else
        @warn "Only-file specified but not found: $(only_file)"
    end
    return only_instances
end

"""
    mark_instance_completed(log_file::String, order::String, scenario::Int, p::Int)

Append a completed instance record to the progress log.
Creates log file with appropriate header if not already present.

# Arguments
- `log_file::String`: Path to progress log file
- `order::String`: Order identifier, e.g., "O1"
- `scenario::Int`: Scenario number
- `p::Int`: Number of parallel machines

# Side Effects
- Creates `log_file` with header comment if it doesn't exist
- Appends new line: `order,scenario,p`
"""
function mark_instance_completed(log_file::String, order::String, scenario::Int, p::Int)
    # Add header if file is new
    if !isfile(log_file)
        open(log_file, "w") do io
            println(io, "# Batch comparison progress log")
            println(io, "# Format: Order,Scenario,P")
            println(io, "# Each line represents a completed instance")
            println(io, "# To restart from scratch, delete this file")
        end
    end
    
    open(log_file, "a") do io
        println(io, "$(order),$(scenario),$(p)")
    end
end

# Main execution
println("=" ^ 80)
println("Comprehensive Batch Comparison: MILP, Split-Solve-Merge, SA, GRASP")
println("=" ^ 80)
println("Configuration:")
println("  SA runs per instance: $(NUM_SA_RUNS)")
println("  Output CSV: $(OUTPUT_CSV)")
println("  Progress log: $(PROGRESS_LOG)")
println("  Error log: $(ERROR_LOG)")
if DRY_RUN_LIMIT !== nothing
    println("  Instance limit: $(DRY_RUN_LIMIT) (dry run mode)")
end
if BETA_OVERRIDE !== nothing
    println("  Beta override: $(BETA_OVERRIDE) (will override values from settings files)")
end
if ONLY_FILE !== nothing
    println("  Only-file: $(ONLY_FILE) (will run only instances listed in file)")
end
if !isempty(SKIP_METHODS)
    println("  Skipped methods: $(join(sort(collect(SKIP_METHODS)), ", "))")
end
println("=" ^ 80)

# Ensure output dir exists
mkpath(OUTPUT_DIR)

# Initialize error log with header
if !isfile(ERROR_LOG)
    open(ERROR_LOG, "w") do io
        timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
        println(io, "# Batch Comparison Error Log")
        println(io, "# Created: ", timestamp)
        println(io, "# This file contains detailed error information for debugging")
        println(io, "# Each error entry includes timestamp, method, instance, error type, message, and stacktrace")
        println(io, "")
    end
else
    # Add session separator if file exists
    open(ERROR_LOG, "a") do io
        timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
        println(io, "\n\n")
        println(io, "#" ^ 80)
        println(io, "# New Session: ", timestamp)
        println(io, "#" ^ 80)
        println(io, "")
    end
end

# Load completed instances from progress log
completed_instances = load_completed_instances(PROGRESS_LOG)

# Load only-instances if specified
only_instances = ONLY_FILE !== nothing ? load_only_instances(ONLY_FILE) : nothing

# Get all E_* and H_* files and group them by (Order, Scenario, P)
settings_dir = datadir("settings")
e_files = filter(f -> startswith(f, "E_") && endswith(f, ".jl"), readdir(settings_dir))
h_files = filter(f -> startswith(f, "H_") && endswith(f, ".jl"), readdir(settings_dir))
sort!(e_files)
sort!(h_files)

println("\nFound $(length(e_files)) E_* instances (exact/MILP)")
println("Found $(length(h_files)) H_* instances (heuristics)")
println()

# Ensure output dir exists
mkpath(OUTPUT_DIR)

# No in-memory results accumulation - write directly to CSV
# Count processed instances for summary
processed_count = 0
method_counts = Dict("MILP" => 0, "SSM" => 0, "SA" => 0, "GRASP" => 0)

# helper: parse and build map
instances = Dict{Tuple{String,Int,Int}, Dict{String,String}}()
for f in e_files
    parsed = parse_instance_name(f)
    if parsed === nothing
        @warn "Could not parse filename: $(f)"
        continue
    end
    key = (parsed.order, parsed.scenario, parsed.p)
    instances[key] = get(instances, key, Dict{String,String}())
    instances[key]["E"] = f
end
for f in h_files
    parsed = parse_instance_name(f)
    if parsed === nothing
        @warn "Could not parse filename: $(f)"
        continue
    end
    key = (parsed.order, parsed.scenario, parsed.p)
    instances[key] = get(instances, key, Dict{String,String}())
    instances[key]["H"] = f
end

"""
    order_index(order::String)::Int

Extract numeric order index from order identifier string.
Used for sorting instances by order number.

# Arguments
- `order::String`: Order identifier, e.g., "O1", "O2", "O10"

# Returns
- `Int`: The numeric value extracted (e.g., 1, 2, 10)
- `0`: If order string does not match pattern "O<number>"

# Examples
```julia
order_index("O1")   # 1
order_index("O10")  # 10
order_index("invalid")  # 0
```
"""
function order_index(order::String)
    m = match(r"O(\d+)", order)
    return m !== nothing ? parse(Int, m.captures[1]) : 0
end

# sort keys by order number, scenario, p
keys_sorted = sort(collect(keys(instances)), by = k -> (order_index(k[1]), k[2], k[3]))

# Apply dry-run limit if set
if DRY_RUN_LIMIT !== nothing
    println("⚠ DRY_RUN_LIMIT set to $(DRY_RUN_LIMIT) — processing only first $(DRY_RUN_LIMIT) instance(s)")
    keys_sorted = keys_sorted[1:min(DRY_RUN_LIMIT, length(keys_sorted))]
end

# Filter by only-file if specified
if only_instances !== nothing
    keys_sorted = filter(k -> k in only_instances, keys_sorted)
    println("⚠ Filtered to $(length(keys_sorted)) instance(s) from only-file")
end

# incremental CSV writer helper
"""
    write_row_incremental(csv_path::String, row_df::DataFrame)

Append a single row (DataFrame) to an existing CSV file or create it if missing.

# Arguments
- `csv_path::String`: Path to CSV file
- `row_df::DataFrame`: Single-row DataFrame to append

# Side Effects
- Creates CSV file if it doesn't exist (overwrites header)
- Appends row to existing CSV file
- Logs warning if write fails

# Behavior
- First call creates file with header from `row_df` column names
- Subsequent calls append row to end of file
"""
function write_row_incremental(csv_path::String, row_df::DataFrame)
    try
        if isfile(csv_path)
            CSV.write(csv_path, row_df; append=true)
        else
            CSV.write(csv_path, row_df)
        end
    catch e
        @warn "Failed to write CSV $(csv_path): $e"
    end
end

total = length(keys_sorted)
global skipped_count = 0
global failed_instances = []

for (idx, key) in enumerate(keys_sorted)
    order, scenario, p = key
    files = instances[key]
    
    # Check if this instance was already completed
    if key in completed_instances
        global skipped_count += 1
        println("⏭  Skipping instance $(idx)/$(total): Order=$(order), Scenario=#$(scenario), p=$(p) [already completed]")
        continue
    end
    
    println("=" ^ 80)
    println("Processing instance $(idx)/$(total): Order=$(order), Scenario=#$(scenario), p=$(p)")
    println("  Files present: ", collect(keys(files)))
    println("=" ^ 80)

    # Wrap entire instance processing in try-catch for maximum robustness
    try
        # Initialize results for this instance
        milp_result = nothing
        sa_result = nothing
        grasp_result = nothing
        ssm_result = nothing

    # If E file exists, run MILP, SA, GRASP
    if haskey(files, "E")
        filename = files["E"]
        filepath = joinpath(settings_dir, filename)
        
        # Try to load the E_ file
        try
            include(filepath)
        catch e
            instance_id = "$(order)_#$(scenario)_$(p)p"
            @warn "Failed to load E_ file $(filename): $e"
            log_error(ERROR_LOG, "E_ File Load", instance_id, e, catch_backtrace())
            println("  ⚠ Skipping all methods for E_ file due to load error")
        end
        
        # Only proceed if order_dict was successfully loaded
        if @isdefined(order_dict)
            # Apply beta override if specified
            if BETA_OVERRIDE !== nothing
                order_dict[:β] = BETA_OVERRIDE
            end
            
            # Run MILP with individual error handling
            if "MILP" in SKIP_METHODS
                println("  Skipping MILP (E_ file)...")
                milp_result = (cost="Skipped", m="Skipped", time="Skipped", error=false)
            else
                println("  Running MILP (E_ file)...")
                try
                    milp_result = run_milp(order_dict)
                    if milp_result.error
                        println("    MILP: ❌ ERROR")
                    else
                        println("    MILP: cost=$(milp_result.cost), m=$(milp_result.m), time=$(round(milp_result.time, digits=3))s")
                    end
                catch e
                    instance_id = "$(order)_#$(scenario)_$(p)p"
                    @warn "Unexpected error running MILP: $e"
                    log_error(ERROR_LOG, "MILP (outer catch)", instance_id, e, catch_backtrace())
                    milp_result = (cost="Err", m="Err", time="Err", error=true)
                    println("    MILP: ❌ ERROR")
                end
            end

            # Run SA with individual error handling
            if "SA" in SKIP_METHODS
                println("  Skipping SA...")
                sa_result = (cost="Skipped", m="Skipped", time="Skipped", error=false)
            else
                println("  Running SA ($(NUM_SA_RUNS) runs)...")
                try
                    sa_result = run_sa_multiple(order_dict, NUM_SA_RUNS)
                    if sa_result.error
                        println("    SA: ❌ ERROR")
                    else
                        println("    SA: cost=$(sa_result.cost), m=$(sa_result.m), time=$(round(sa_result.time, digits=3))s")
                    end
                catch e
                    instance_id = "$(order)_#$(scenario)_$(p)p"
                    @warn "Unexpected error running SA: $e"
                    log_error(ERROR_LOG, "SA (outer catch)", instance_id, e, catch_backtrace())
                    sa_result = (cost="Err", m="Err", time="Err", error=true)
                    println("    SA: ❌ ERROR")
                end
            end

            # Run GRASP with individual error handling
            if "GRASP" in SKIP_METHODS
                println("  Skipping GRASP...")
                grasp_result = (cost="Skipped", m="Skipped", time="Skipped", error=false)
            else
                println("  Running GRASP...")
                try
                    grasp_result = run_grasp_once(order_dict)
                    if grasp_result.error
                        println("    GRASP: ❌ ERROR")
                    else
                        println("    GRASP: cost=$(grasp_result.cost), m=$(grasp_result.m), time=$(round(grasp_result.time, digits=3))s")
                    end
                catch e
                    instance_id = "$(order)_#$(scenario)_$(p)p"
                    @warn "Unexpected error running GRASP: $e"
                    log_error(ERROR_LOG, "GRASP (outer catch)", instance_id, e, catch_backtrace())
                    grasp_result = (cost="Err", m="Err", time="Err", error=true)
                    println("    GRASP: ❌ ERROR")
                end
            end
        end
    end

    # If H file exists, run Split-Solve-Merge
    if haskey(files, "H")
        filename = files["H"]
        filepath = joinpath(settings_dir, filename)
        
        # Try to load the H_ file
        try
            include(filepath)
        catch e
            instance_id = "$(order)_#$(scenario)_$(p)p"
            @warn "Failed to load H_ file $(filename): $e"
            log_error(ERROR_LOG, "H_ File Load", instance_id, e, catch_backtrace())
            println("  ⚠ Skipping Split-Solve-Merge due to load error")
        end
        
        # Only proceed if order_dict was successfully loaded
        if @isdefined(order_dict)
            # Apply beta override if specified
            if BETA_OVERRIDE !== nothing
                order_dict[:β] = BETA_OVERRIDE
            end
            
            if "SSM" in SKIP_METHODS
                println("  Skipping Split-Solve-Merge (H_ file)...")
                ssm_result = (cost="Skipped", m="Skipped", time="Skipped", error=false)
            else
                println("  Running Split-Solve-Merge (H_ file)...")
                try
                    ssm_result = run_split_solve_merge(order_dict)
                    if ssm_result.error
                        println("    SSM: ❌ ERROR")
                    else
                        println("    SSM: cost=$(ssm_result.cost), m=$(ssm_result.m), time=$(round(ssm_result.time, digits=3))s")
                    end
                catch e
                    instance_id = "$(order)_#$(scenario)_$(p)p"
                    @warn "Unexpected error running Split-Solve-Merge: $e"
                    log_error(ERROR_LOG, "SSM (outer catch)", instance_id, e, catch_backtrace())
                    ssm_result = (cost="Err", m="Err", time="Err", error=true)
                    println("    SSM: ❌ ERROR")
                end
            end
        end
    end

    # Combine results into single row for this instance
    row = (
        Order = order,
        Scenario = scenario,
        P = p,
        MILP_Cost = milp_result !== nothing ? milp_result.cost : missing,
        MILP_Slots = milp_result !== nothing ? milp_result.m : missing,
        MILP_Time = milp_result !== nothing ? milp_result.time : missing,
        SSM_Cost = ssm_result !== nothing ? ssm_result.cost : missing,
        SSM_Slots = ssm_result !== nothing ? ssm_result.m : missing,
        SSM_Time = ssm_result !== nothing ? ssm_result.time : missing,
        SA_Cost = sa_result !== nothing ? sa_result.cost : missing,
        SA_Slots = sa_result !== nothing ? sa_result.m : missing,
        SA_Time = sa_result !== nothing ? sa_result.time : missing,
        GRASP_Cost = grasp_result !== nothing ? grasp_result.cost : missing,
        GRASP_Slots = grasp_result !== nothing ? grasp_result.m : missing,
        GRASP_Time = grasp_result !== nothing ? grasp_result.time : missing
    )
    # Write to CSV incrementally
    df_row = DataFrame([row])
    write_row_incremental(OUTPUT_CSV, df_row)
    # also write to per-scenario raw CSV for safety
    scenario_raw = joinpath(OUTPUT_DIR, "raw_comparison_scenario_$(scenario).csv")
    write_row_incremental(scenario_raw, df_row)
    
    # Update counters
    global processed_count += 1
    !ismissing(row.MILP_Cost) && row.MILP_Cost != "Skipped" && (method_counts["MILP"] += 1)
    !ismissing(row.SSM_Cost) && row.SSM_Cost != "Skipped" && (method_counts["SSM"] += 1)
    !ismissing(row.SA_Cost) && row.SA_Cost != "Skipped" && (method_counts["SA"] += 1)
    !ismissing(row.GRASP_Cost) && row.GRASP_Cost != "Skipped" && (method_counts["GRASP"] += 1)        # Mark instance as completed in progress log
        mark_instance_completed(PROGRESS_LOG, order, scenario, p)
        
        # MEMORY MANAGEMENT: Force garbage collection every 10 instances
        if idx % 10 == 0
            GC.gc()
            mem_mb = Base.gc_live_bytes() / 1024^2
            println("  [Memory: $(round(mem_mb, digits=1)) MB live after $(idx)/$(total) instances]")
        end

    catch e
        # Catastrophic failure for this instance - log and continue
        instance_id = "$(order)_#$(scenario)_$(p)p"
        @error "Catastrophic failure processing instance $(instance_id): $e"
        log_error(ERROR_LOG, "CATASTROPHIC INSTANCE FAILURE", instance_id, e, catch_backtrace())
        println("  ❌ INSTANCE FAILED - continuing with next instance")
        global failed_instances
        push!(failed_instances, (order=order, scenario=scenario, p=p, error=string(e)))
        
        # Write a row with all errors to maintain CSV structure
        try
            row = (
                Order = order,
                Scenario = scenario,
                P = p,
                MILP_Cost = "Err",
                MILP_Slots = "Err",
                MILP_Time = "Err",
                SSM_Cost = "Err",
                SSM_Slots = "Err",
                SSM_Time = "Err",
                SA_Cost = "Err",
                SA_Slots = "Err",
                SA_Time = "Err",
                GRASP_Cost = "Err",
                GRASP_Slots = "Err",
                GRASP_Time = "Err"
            )
            df_row = DataFrame([row])
            write_row_incremental(OUTPUT_CSV, df_row)
        catch csv_error
            @warn "Could not write error row to CSV: $(csv_error)"
        end
    end

    println()
end

# Print skip summary
if skipped_count > 0
    println("=" ^ 80)
    println("⏭  Skipped $(skipped_count) already completed instance(s)")
    println("=" ^ 80)
end

# Print failed instances summary
if !isempty(failed_instances)
    println("\n" * "=" ^ 80)
    println("⚠  FAILED INSTANCES SUMMARY")
    println("=" ^ 80)
    println("Total failed instances: $(length(failed_instances))")
    for (i, failure) in enumerate(failed_instances)
        println("  $(i). $(failure.order)_#$(failure.scenario)_$(failure.p)p")
        println("     Error: $(failure.error)")
    end
    println("=" ^ 80)
end

# Print summary statistics
println("\n" * "=" ^ 80)
println("SUMMARY STATISTICS")
println("=" ^ 80)
println(@sprintf("Total instances processed: %d", processed_count))
println(@sprintf("Instances with MILP results: %d", method_counts["MILP"]))
println(@sprintf("Instances with SSM results: %d", method_counts["SSM"]))
println(@sprintf("Instances with SA results: %d", method_counts["SA"]))
println(@sprintf("Instances with GRASP results: %d", method_counts["GRASP"]))

# Final memory report
mem_mb = Base.gc_live_bytes() / 1024^2
println(@sprintf("Final memory usage: %.1f MB", mem_mb))

# Check if error log has content beyond header
error_log_has_errors = false
if isfile(ERROR_LOG)
    error_log_size = filesize(ERROR_LOG)
    if error_log_size > 500  # More than just header
        error_log_has_errors = true
    end
end

if error_log_has_errors || !isempty(failed_instances)
    println("\n⚠  ERRORS OCCURRED - detailed error information written to:")
    println("   $(ERROR_LOG)")
end

println("\nNote: For detailed method comparisons, run generate_latex_tables.jl")

println("=" ^ 80)
