module Grasp
export run, create_grasp_partition, run_greedy

using DrWatson
@quickactivate "ShoeOptSetupTime"

using Logging, Dates, DataFrames, Random, Printf

# Import needed functions from SimulatedAnnealing
include(scriptsdir("simulated_annealing.jl"))
using .SimulatedAnnealing: split_quantity_randomly, get_cost_from_shelves, 
                                partition_to_dataframe, extract_slots_from_partition

"""
    even_split(total::Int, parts::Int)::Vector{Int}

Deterministically split a total quantity into parts with equal distribution.

# Arguments
- `total::Int`: Total quantity to be divided
- `parts::Int`: Number of parts to split into (must be positive)

# Returns
- `Vector{Int}`: Vector of `parts` integers that sum to `total`, distributed as evenly as possible

# Algorithm
Uses floor division with remainder distribution. When `total % parts != 0`, earlier parts 
receive one extra unit each to account for the remainder. This ensures deterministic behavior.

# Examples
```julia
even_split(10, 3)  # [4, 3, 3]
even_split(7, 2)   # [4, 3]
even_split(5, 1)   # [5]
```

# Errors
Throws an error if `parts <= 0`
"""
function even_split(total::Int, parts::Int)
    parts <= 0 && error("parts must be positive")
    parts == 1 && return [total]
    base = total ÷ parts
    remn = total % parts
    v = fill(base, parts)
    for i in 1:remn
        v[i] += 1
    end
    return v
end

"""
    create_all_subjobs(g, n, o; split_method::Symbol = :random)::Vector

Create the complete list of subjobs for all jobs in the scheduling problem.

Each job may require multiple molds, resulting in multiple subjobs. This function 
decomposes jobs into their constituent subjobs for granular scheduling.

# Arguments
- `g`: Vector of job identifiers
- `n`: Vector of total quantities for each job (parallel to `g`)
- `o`: Vector indicating number of molds required per job (parallel to `g`)
- `split_method::Symbol`: Method for dividing job quantities across molds
  - `:random` - Uses stochastic split (SimulatedAnnealing.split_quantity_randomly)
  - `:even` - Uses deterministic equal distribution (even_split)

# Returns
- `Vector`: List of NamedTuples with fields `(jobid, moldid, quantity)`

# Behavior
- Single-mold jobs (o[i] == 1): Create one subjob with full quantity
- Multi-mold jobs (o[i] > 1): Split quantity across molds using specified method

# Examples
```julia
g = [1, 2]
n = [10, 5]
o = [2, 1]
create_all_subjobs(g, n, o; split_method=:even)
# Returns subjobs like: [(jobid=1, moldid=1, quantity=5), 
#                        (jobid=1, moldid=2, quantity=5),
#                        (jobid=2, moldid=1, quantity=5)]
```
"""
function create_all_subjobs(g, n, o; split_method::Symbol = :random)
    all_subjobs = []
    for (idx, jobid) in enumerate(g)
        total_qty = n[idx]
        num_molds = o[idx]
        if num_molds == 1
            push!(all_subjobs, (jobid=jobid, moldid=1, quantity=total_qty))
        else
            mold_quantities = split_method == :random ?
                split_quantity_randomly(total_qty, num_molds) :
                even_split(total_qty, num_molds)
            for (moldid, qty) in enumerate(mold_quantities)
                push!(all_subjobs, (jobid=jobid, moldid=moldid, quantity=qty))
            end
        end
    end
    return all_subjobs
end

"""
    select_subjob(subjobs; mode::Symbol = :probabilistic)::Tuple

Select a subjob from the candidate list using specified selection strategy.

# Arguments
- `subjobs`: Vector of subjob NamedTuples (each with fields: jobid, moldid, quantity)
- `mode::Symbol`: Selection strategy
  - `:probabilistic` - Roulette-wheel selection proportional to quantity (GRASP behavior)
  - `:greedy` - Deterministically select subjob with maximum quantity (ties: first occurrence)

# Returns
- `Tuple{NamedTuple, Vector}`: (selected_subjob, remaining_subjobs)

# Selection Strategies

## Probabilistic Mode
Uses roulette-wheel selection where selection probability is proportional to subjob quantity:
- P(subjob_i) = quantity_i / sum(all quantities)
- Introduces randomness for GRASP exploration

## Greedy Mode
Deterministic selection of the largest subjob:
- Always selects subjob with maximum quantity
- Ties resolved by first occurrence in the list
- Used for pure greedy construction

# Errors
Throws error if `subjobs` is empty

# Examples
```julia
subjobs = [(jobid=1, moldid=1, quantity=10), 
           (jobid=2, moldid=1, quantity=5)]
selected, remaining = select_subjob(subjobs; mode=:greedy)
# selected = (jobid=1, moldid=1, quantity=10)
# remaining = [(jobid=2, moldid=1, quantity=5)]
```
"""
function select_subjob(subjobs; mode::Symbol = :probabilistic)
    if isempty(subjobs)
        error("select_subjob called with empty list")
    end
    if mode == :greedy
        # Deterministic: pick maximum quantity (first occurrence on ties)
        quantities = [sj.quantity for sj in subjobs]
        idx = findfirst(==(maximum(quantities)), quantities)
        selected = subjobs[idx]
        remaining = [subjobs[i] for i in 1:length(subjobs) if i != idx]
        return selected, remaining
    end

    total_qty = sum(sj.quantity for sj in subjobs)
    if total_qty == 0
        idx = rand(1:length(subjobs))
        selected = subjobs[idx]
        remaining = [subjobs[i] for i in 1:length(subjobs) if i != idx]
        return selected, remaining
    end
    probabilities = [sj.quantity / total_qty for sj in subjobs]
    cumsum_probs = cumsum(probabilities)
    r = rand()
    idx = findfirst(>=(r), cumsum_probs)
    selected = subjobs[idx]
    remaining = [subjobs[i] for i in 1:length(subjobs) if i != idx]
    return selected, remaining
end

"""
    get_shelf_loads(shelves)::Vector{Int}

Calculate the total load (sum of quantities) for each shelf.

# Arguments
- `shelves`: Vector of shelf vectors, where each shelf contains subjob NamedTuples

# Returns
- `Vector{Int}`: Load for each shelf (sum of all subjob quantities on that shelf)

# Note
Empty shelves have a load of 0.
"""
function get_shelf_loads(shelves)
    loads = zeros(Int, length(shelves))
    for (shelf_idx, shelf) in enumerate(shelves)
        if !isempty(shelf)
            loads[shelf_idx] = sum(sj.quantity for sj in shelf)
        end
    end
    return loads
end

"""
    place_subjob_greedy(subjob, shelves)::Vector

Place a subjob on the shelf that results in minimum maximum load (load-balancing heuristic).

# Arguments
- `subjob`: NamedTuple with fields (jobid, moldid, quantity)
- `shelves`: Current shelf configuration (vector of vectors of subjobs)

# Returns
- `Vector`: New shelf configuration with subjob assigned to optimal shelf

# Algorithm
Implements a greedy load-balancing strategy:
1. Calculate current load of each shelf
2. Compute hypothetical load if subjob is added to each shelf
3. Assign subjob to shelf that minimizes the maximum load
4. Return deep copy of updated shelf configuration

# Note
This is the core constructive heuristic for GRASP - it greedily balances load 
across parallel machines while building the schedule.
"""
function place_subjob_greedy(subjob, shelves)
    new_shelves = deepcopy(shelves)
    loads = get_shelf_loads(new_shelves)
    new_loads = loads .+ subjob.quantity
    min_load_shelf = argmin(new_loads)
    push!(new_shelves[min_load_shelf], subjob)
    return new_shelves
end

"""
    create_grasp_partition(g, n, o, p; selection::Symbol = :probabilistic, split_method::Symbol = :random)::Vector

Construct a single partition (schedule) using GRASP methodology.

# Arguments
- `g`: Vector of job identifiers
- `n`: Vector of job quantities (parallel to `g`)
- `o`: Vector of molds per job (parallel to `g`)
- `p`: Number of parallel shelves (machines)
- `selection::Symbol`: Subjob selection strategy
  - `:probabilistic` - Roulette-wheel selection (randomized, for GRASP)
  - `:greedy` - Deterministic maximum-quantity selection
- `split_method::Symbol`: Quantity splitting method for multi-mold jobs
  - `:random` - Stochastic split
  - `:even` - Deterministic equal distribution

# Returns
- `Vector{Vector}`: Partition represented as `p` shelves, each containing assigned subjobs

# Algorithm
1. Initialize `p` empty shelves
2. Create all subjobs by splitting jobs across molds
3. Iteratively select subjobs using specified strategy
4. Place each selected subjob on the shelf that minimizes load imbalance
5. Continue until all subjobs are assigned

# Randomness vs Determinism
- **Pure GRASP**: `selection=:probabilistic`, `split_method=:random` (default)
- **Pure Greedy**: `selection=:greedy`, `split_method=:even` (deterministic)

# See Also
- [`select_subjob`](@ref) for selection strategies
- [`place_subjob_greedy`](@ref) for placement heuristic
"""
function create_grasp_partition(g, n, o, p; selection::Symbol = :probabilistic, split_method::Symbol = :random)
    shelves = [[] for _ in 1:p]
    remaining_subjobs = create_all_subjobs(g, n, o; split_method=split_method)
    while !isempty(remaining_subjobs)
        selected_subjob, remaining_subjobs = select_subjob(remaining_subjobs; mode=selection)
        shelves = place_subjob_greedy(selected_subjob, shelves)
    end
    return shelves
end

"""
    run(order_dict; log_every=10, selection=:probabilistic, split_method=:random, iterations=nothing)::Vector

Execute GRASP (Greedy Randomized Adaptive Search Procedure) for parallel machine scheduling.

# Arguments
- `order_dict::Dict`: Problem instance dictionary with required keys:
  - `:g` - Job identifiers
  - `:n` - Job quantities
  - `:o` - Molds per job
  - `:p` - Number of parallel shelves/machines
  - `:α` - Shelf penalty coefficient
  - `:β` - Job delay penalty coefficient
  - `:Nit` - Default number of iterations (overridden by `iterations` argument)
- `log_every::Int=10`: Logging frequency (iterations per log message, 0 to disable)
- `selection::Symbol=:probabilistic`: Subjob selection strategy (`:probabilistic` or `:greedy`)
- `split_method::Symbol=:random`: Mold quantity splitting method (`:random` or `:even`)
- `iterations::Union{Nothing,Int}=nothing`: Override iteration count (uses `order_dict[:Nit]` if `nothing`)

# Returns
- `Vector{Vector}`: Best partition found (shelves with assigned subjobs)

# Algorithm
Multi-start GRASP procedure:
1. For each iteration:
   - Construct a new partition using randomized greedy heuristic
   - Evaluate cost using objective function
   - Track best solution found
2. Return partition with minimum cost

# Output
Logs progress at specified intervals and displays final results including:
- Best cost, shelf count, and maximum load
- Detailed partition structure
- DataFrames with slot assignments

# GRASP vs Greedy Mode
- **GRASP** (default): `selection=:probabilistic`, `split_method=:random` → Multiple diverse solutions
- **Pure Greedy**: `selection=:greedy`, `split_method=:even`, `iterations=1` → Single deterministic solution

# Examples
```julia
# Standard GRASP with 100 iterations
shelves = run(order_dict)

# Pure greedy (deterministic)
shelves = run(order_dict; selection=:greedy, split_method=:even, iterations=1)

# Quick test with 5 iterations
shelves = run(order_dict; iterations=5, log_every=1)
```

# See Also
- [`run_greedy`](@ref) for dedicated greedy solver
- [`create_grasp_partition`](@ref) for single partition construction
"""
function run(order_dict; log_every::Int=10, selection::Symbol=:probabilistic, split_method::Symbol=:random, iterations::Union{Nothing,Int}=nothing)
    @info "Running GRASP with order_dict:"
    @info order_dict
    @unpack g, n, o, p, α, β, Nit = order_dict
    iters = isnothing(iterations) ? Nit : iterations
    best_shelves = nothing
    best_max_sum = Inf
    best_m = Inf
    best_cost = Inf
    @info @sprintf("GRASP params: iterations=%d selection=%s split=%s", iters, String(Symbol(selection)), String(Symbol(split_method)))
    for iteration in 1:iters
        current_shelves = create_grasp_partition(g, n, o, p; selection=selection, split_method=split_method)
        current_max_sum, current_m, current_cost = get_cost_from_shelves(current_shelves, α, β)
        if current_cost < best_cost
            best_cost = current_cost
            best_max_sum = current_max_sum
            best_m = current_m
            best_shelves = deepcopy(current_shelves)
        end
        if log_every > 0 && iteration % log_every == 0
            @info @sprintf(
                "it=%5d curr_cost=%.6f best_cost=%.6f",
                iteration, current_cost, best_cost)
        end
    end
    @info @sprintf("GRASP finished: iterations=%d best_cost=%.6f best_m=%d best_max_sum=%.6f", 
                   iters, best_cost, best_m, best_max_sum)
    @info "Best partition (shelves) snapshot:"
    for (i, shelf) in enumerate(best_shelves)
        println("  Shelf $i: ", join(["(job=$(j.jobid), mold=$(j.moldid), qty=$(j.quantity))" for j in shelf], ", "))
    end
    df = partition_to_dataframe(best_shelves)
    @info df
    dfp = extract_slots_from_partition(df)
    @info dfp
    return best_shelves
end

"""
    run_greedy(order_dict; log_every=0)::NamedTuple

Compute a single deterministic greedy solution for single-mold scenarios.

# Arguments
- `order_dict::Dict`: Problem instance dictionary (see [`run`](@ref) for required keys)
- `log_every::Int=0`: Logging level (≥ 0 logs summary, < 0 silent)

# Returns
- `NamedTuple` with fields:
  - `shelves`: Optimal greedy partition
  - `max_sum`: Maximum shelf load
  - `m`: Number of non-empty shelves
  - `cost`: Total objective cost
  - `time`: Execution time in seconds

# Restrictions
**Single-mold jobs only**: All jobs must have `o[j] == 1` (one mold per job).
Throws error if multi-mold jobs are present.

# Algorithm
Uses deterministic greedy construction:
- `selection=:greedy` - Always select largest remaining subjob
- `split_method=:even` - Deterministic quantity splitting (though not needed for single-mold)
- Single iteration (no randomness)

# Use Cases
- Quick baseline solution for comparison
- Deterministic lower bound for metaheuristics
- Single-mold problem instances only

# Errors
Throws error if any job requires multiple molds (`o[j] != 1`)

# Examples
```julia
result = run_greedy(order_dict)
# Returns: (shelves=[...], max_sum=150, m=3, cost=450.0, time=0.002)
```

# See Also
- [`run`](@ref) for full GRASP with multi-mold support
"""
function run_greedy(order_dict; log_every::Int=0)
    @unpack g, n, o, p, α, β = order_dict
    # Enforce single-mold scenarios: all o[j] must be 1
    if any(x -> x != 1, o)
        error("Greedy only available for single-mold scenarios (all o[j] == 1): EXIT")
    end
    elapsed = @elapsed begin
        shelves = create_grasp_partition(g, n, o, p; selection=:greedy, split_method=:even)
        global _rg_shelves = shelves
    end
    shelves = _rg_shelves
    max_sum, m, cost = get_cost_from_shelves(shelves, α, β)
    if log_every >= 0
        @info @sprintf("Greedy finished: cost=%.6f m=%d max_sum=%.6f time=%.3f", cost, m, max_sum, elapsed)
    end
    return (shelves=shelves, max_sum=max_sum, m=m, cost=cost, time=elapsed)
end

end
