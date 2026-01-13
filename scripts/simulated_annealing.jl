#=  Copyright (C) 2024
Nuno David Lopes.
Created:  2024/04/09
Last changed - N. Lopes: 2026/01/13 11:57:46
=#

module SimulatedAnnealing
export run, create_random_partition, partition_to_dataframe
export extract_slots_from_partition, neighbor_ppartition, neighbor_quantity_redistribution
export get_cost_from_shelves, get_cost_from_slots, split_quantity_randomly, two_step_neighbor

using DrWatson
@quickactivate "ShoeOptSetupTime"

using Logging, Dates, DataFrames, Random, Combinatorics, Printf

"""
    split_quantity_randomly(total_quantity::Int, num_parts::Int; method::Symbol=:stick_breaking)::Vector{Int}

Randomly split a total quantity into a specified number of non-negative integer parts.

# Arguments
- `total_quantity::Int`: Total quantity to be divided (must be non-negative)
- `num_parts::Int`: Number of parts to split into (must be positive)
- `method::Symbol=:stick_breaking`: Splitting method
  - `:stick_breaking` - Fast uniform distribution over compositions (default)
  - `:uniform` - Exact uniform selection from all compositions (slower for large inputs)

# Returns
- `Vector{Int}`: Vector of `num_parts` non-negative integers summing to `total_quantity`

# Methods

## Stick-Breaking Method (default)
Uses the "stars and bars" technique with random permutation. Computationally efficient
and produces a uniform distribution over all integer compositions.

## Uniform Method
Enumerates all possible compositions and selects uniformly. Exact but computationally
expensive for large `total_quantity` and `num_parts`.

# Examples
```julia
split_quantity_randomly(10, 3)  # e.g., [3, 5, 2]
split_quantity_randomly(10, 3; method=:uniform)  # Exact uniform sampling
split_quantity_randomly(0, 2)   # [0, 0]
split_quantity_randomly(5, 1)   # [5]
```

# Errors
- Throws error if `num_parts <= 0`
- Throws error if `total_quantity < 0`
- Throws error for unknown `method`
"""
function split_quantity_randomly(total_quantity::Int, num_parts::Int; 
                                  method::Symbol=:stick_breaking)
    if num_parts <= 0
        error("num_parts must be positive")
    end
    if total_quantity < 0
        error("total_quantity must be non-negative")
    end
    if num_parts == 1
        return [total_quantity]
    end
    if total_quantity == 0
        return zeros(Int, num_parts)
    end
    if method == :stick_breaking
        return _split_stick_breaking(total_quantity, num_parts)
    elseif method == :uniform
        return _split_uniform(total_quantity, num_parts)
    else
        error("Unknown method: $method. Use :stick_breaking or :uniform")
    end
end

"""
    _split_stick_breaking(total_quantity::Int, num_parts::Int)::Vector{Int}

Internal function implementing the stick-breaking (stars and bars) method for random composition.

Uses random permutation to select divider positions, ensuring uniform distribution over
all integer compositions.

# Algorithm
1. Generate `total_quantity + num_parts` positions
2. Randomly select `num_parts - 1` divider positions
3. Compute part sizes from gaps between dividers
4. Subtract 1 from each part to account for offset

# Note
This is an internal helper function. Use [`split_quantity_randomly`](@ref) instead.
"""
function _split_stick_breaking(total_quantity::Int, num_parts::Int)
    n_plus_k = total_quantity + num_parts
    all_positions = randperm(n_plus_k - 1)[1:(num_parts - 1)]
    cut_points = sort(all_positions)
    quantities = Vector{Int}(undef, num_parts)
    prev = 0
    for i in 1:num_parts-1
        quantities[i] = cut_points[i] - prev
        prev = cut_points[i]
    end
    quantities[end] = n_plus_k - prev
    return quantities .- 1
end

"""
    _split_uniform(total_quantity::Int, num_parts::Int)::Vector{Int}

Internal function implementing exact uniform sampling from all integer compositions.

Enumerates all possible ways to partition `total_quantity` into `num_parts` non-negative
integers and selects one uniformly at random.

# Warning
Computational complexity grows rapidly with input size. Use `:stick_breaking` method
for large inputs.

# Note
This is an internal helper function. Use [`split_quantity_randomly`](@ref) instead.
"""
function _split_uniform(total_quantity::Int, num_parts::Int)
    n_plus_k = total_quantity + num_parts
    all_compositions = Vector{Vector{Int}}()
    for dividers in combinations(1:(n_plus_k-1), num_parts-1)
        comp = zeros(Int, num_parts)
        prev = 0
        for i in 1:num_parts-1
            comp[i] = dividers[i] - prev
            prev = dividers[i]
        end
        comp[num_parts] = n_plus_k - prev
        push!(all_compositions, comp .- 1)
    end
    return rand(all_compositions)
end

"""
    create_random_partition(g, n, o, p)::Vector{Vector}

Generate a random initial partition of jobs across parallel shelves.

# Arguments
- `g`: Vector of job identifiers
- `n`: Vector of job quantities (parallel to `g`)
- `o`: Vector of molds per job (parallel to `g`)
- `p`: Number of parallel shelves (machines)

# Returns
- `Vector{Vector}`: Partition represented as `p` shelves, each containing randomly assigned subjobs

# Algorithm
1. For each job:
   - Split job quantity across its molds using random composition
   - Create subjobs (jobid, moldid, quantity) for each mold
   - Randomly assign each subjob to one of the `p` shelves
2. Return shelf configuration

# Structure
Each subjob is a NamedTuple: `(jobid=Int, moldid=Int, quantity=Int)`

# Use Cases
- Initial solution for Simulated Annealing
- Random restart for metaheuristics
- Baseline comparison

# Examples
```julia
g = [1, 2]
n = [10, 5]
o = [2, 1]
p = 3
shelves = create_random_partition(g, n, o, p)
# Returns 3 shelves with randomly assigned subjobs
```
"""
function create_random_partition(g, n, o, p)
    Ps = [shelf for shelf ∈ 1:p]
    shelves = [[] for _ ∈ 1:p]
    for (idx, jobid) ∈ enumerate(g)
        total_qty = n[idx]
        num_molds = o[idx]
        mold_quantities = split_quantity_randomly(total_qty, num_molds)
        for (moldid, qty) ∈ enumerate(mold_quantities)
            shelf = rand(Ps)
            push!(shelves[shelf], (jobid=jobid, moldid=moldid, quantity=qty))
        end
    end
    return shelves
end

"""
    get_cost_from_shelves(shelves, α, β)::Tuple{Int, Int, Float64}

Calculate objective cost from shelf partition using problem-specific cost function.

# Arguments
- `shelves`: Partition as vector of shelves, each containing subjobs
- `α`: Penalty coefficient for maximum shelf load
- `β`: Penalty coefficient for number of distinct shelf levels

# Returns
- `Tuple{Int, Int, Float64}`: (max_sum, m, cost)
  - `max_sum::Int`: Maximum cumulative load across all shelves
  - `m::Int`: Number of distinct shelf levels (unique cumulative quantities)
  - `cost::Float64`: Total objective cost = α × max_sum + β × m

# Cost Function
The objective balances two competing goals:
- Minimize maximum shelf load (`max_sum`) - workload balancing
- Minimize number of shelf levels (`m`) - reduce interruptions

# Algorithm
1. Compute cumulative quantities for all positions in all shelves
2. Find maximum cumulative quantity (`max_sum`)
3. Count unique cumulative values (`m`)
4. Calculate weighted cost

# Edge Cases
- Empty shelves return (0, 0, 0)
- Shelves with no subjobs return (0, 0, 0)

# See Also
- [`get_cost_from_slots`](@ref) for slot-based cost calculation
"""
function get_cost_from_shelves(shelves, α, β)
    isempty(shelves) && return 0, 0
    all_cums = Int[]
    max_cum = 0
    for shelf in shelves
        running = 0
        for job in shelf
            running += job.quantity
            push!(all_cums, running)
        end
        running > max_cum && (max_cum = running)
    end
    if isempty(all_cums)
        return 0, 0
    end
    m = length(unique(all_cums))
    cost = α * max_cum + β * m
    return max_cum, m, cost
end

"""
    partition_to_dataframe(shelves)::DataFrame

Convert shelf partition to a wide-format DataFrame for analysis and visualization.

# Arguments
- `shelves`: Partition as vector of shelves, each containing subjobs

# Returns
- `DataFrame`: Wide format with columns for each shelf:
  - `S<i>JobID`: Job identifier for shelf `i`
  - `S<i>MoldID`: Mold identifier for shelf `i`
  - `S<i>Qty`: Quantity for shelf `i`

# Structure
- Each row represents a position in the shelf sequence
- Zero values indicate empty positions
- Number of columns = 3 × number of shelves

# Use Cases
- Data export for analysis
- Visualization preparation
- Integration with slot extraction

# Examples
```julia
shelves = [[(jobid=1, moldid=1, quantity=5)], 
           [(jobid=2, moldid=1, quantity=3)]]
df = partition_to_dataframe(shelves)
# DataFrame with columns: S1JobID, S1MoldID, S1Qty, S2JobID, S2MoldID, S2Qty
```

# See Also
- [`extract_slots_from_partition`](@ref) for converting to slot representation
"""
function partition_to_dataframe(shelves)
    max_subjobs = sum(length, shelves)
    p = length(shelves)
    jobid_matrix = fill(0, max_subjobs, p)
    moldid_matrix = fill(0, max_subjobs, p)
    quantity_matrix = fill(0, max_subjobs, p)
    for (shelf_idx, shelf_subjobs) ∈ enumerate(shelves)
        for (row_idx, subjob_info) ∈ enumerate(shelf_subjobs)
            jobid_matrix[row_idx, shelf_idx] = subjob_info.jobid
            moldid_matrix[row_idx, shelf_idx] = subjob_info.moldid
            quantity_matrix[row_idx, shelf_idx] = subjob_info.quantity
        end
    end
    df_cols = []
    col_names = String[]
    for i in 1:p
        push!(df_cols, jobid_matrix[:, i])
        push!(col_names, "S$(i)JobID")
        push!(df_cols, moldid_matrix[:, i])
        push!(col_names, "S$(i)MoldID")
        push!(df_cols, quantity_matrix[:, i])
        push!(col_names, "S$(i)Qty")
    end
    dfp = DataFrame(df_cols, col_names)
    return dfp
end

"""
    extract_slots_from_partition(dfp::DataFrame)::DataFrame

Convert partition DataFrame to slot-based representation where each slot has uniform processing.

# Arguments
- `dfp::DataFrame`: Partition in wide format (output of [`partition_to_dataframe`](@ref))

# Returns
- `DataFrame`: Slot representation where each row is a processing slot
  - Each slot has the same structure as input but represents synchronized processing
  - Within each slot, all parallel machines process simultaneously until first completion

# Algorithm
Implements the "simultaneous interruption" mechanism:
1. For each row in partition, find minimum quantity across active shelves
2. Create a slot with this minimum quantity for all active shelves
3. Reduce remaining quantities by the minimum
4. Repeat until all quantities are processed

# Simultaneous Interruption Model
This represents the key problem constraint: when the first job in a parallel group
completes, all other jobs must interrupt, creating a new processing slot.

# Examples
```julia
# Partition: Shelf 1 has job with qty=5, Shelf 2 has job with qty=3
dfp = partition_to_dataframe(shelves)
dfs = extract_slots_from_partition(dfp)
# Result: 2 slots
#   Slot 1: both shelves process qty=3 (min)
#   Slot 2: only shelf 1 processes remaining qty=2
```

# See Also
- [`partition_to_dataframe`](@ref) for creating input format
- [`get_cost_from_slots`](@ref) for cost calculation on slots
"""
function extract_slots_from_partition(dfp)
    nrows, ncols = size(dfp)
    p = ncols ÷ 3
    result_rows = []
    for row in 1:nrows
        working_row = [dfp[row, col] for col in 1:ncols]
        while true
            quantities = []
            for shelf in 1:p
                qty_col = 3 * shelf
                qty = working_row[qty_col]
                if qty > 0
                    push!(quantities, qty)
                end
            end
            if isempty(quantities)
                break
            end
            min_qty = minimum(quantities)
            slot_row = zeros(Int, ncols)
            for shelf in 1:p
                jobid_col = 3 * shelf - 2
                moldid_col = 3 * shelf - 1
                qty_col = 3 * shelf
                qty = working_row[qty_col]
                jid = working_row[jobid_col]
                mid = working_row[moldid_col]
                if qty > 0
                    slot_row[jobid_col] = jid
                    slot_row[moldid_col] = mid
                    slot_row[qty_col] = min_qty
                    working_row[qty_col] = qty - min_qty
                end
            end
            push!(result_rows, slot_row)
        end
    end
    if isempty(result_rows)
        return DataFrame([Int[] for _ in 1:ncols], names(dfp))
    end
    result_matrix = vcat([row' for row in result_rows]...)
    dfs = DataFrame(result_matrix, names(dfp))
    return dfs
end

"""
    get_cost_from_slots(dfs::DataFrame, α, β)::Tuple{Int, Int, Float64}

Calculate objective cost from slot-based representation.

# Arguments
- `dfs::DataFrame`: Slot representation (output of [`extract_slots_from_partition`](@ref))
- `α`: Penalty coefficient for maximum shelf load
- `β`: Penalty coefficient for number of slots

# Returns
- `Tuple{Int, Int, Float64}`: (max_sum, m, cost)
  - `max_sum::Int`: Maximum total load on any shelf
  - `m::Int`: Number of non-empty slots (shelf levels)
  - `cost::Float64`: Total objective cost = α × max_sum + β × m

# Difference from get_cost_from_shelves
This function operates on slot representation where:
- `m` counts processing slots (interruption points)
- `max_sum` is total load per shelf (column sum)

# Edge Cases
- Empty DataFrame returns (0, 0, 0)

# See Also
- [`get_cost_from_shelves`](@ref) for partition-based cost calculation
- [`extract_slots_from_partition`](@ref) for creating slot representation
"""
function get_cost_from_slots(dfs, α, β)
    nrows, ncols = size(dfs)
    if nrows == 0 || ncols == 0
        return 0, 0, 0
    end
    p = ncols ÷ 3
    column_sums = Int[]
    for shelf_idx in 1:p
        qty_col = 3 * shelf_idx
        col_sum = sum(dfs[:, qty_col])
        push!(column_sums, col_sum)
    end
    max_cum = isempty(column_sums) ? 0 : maximum(column_sums)
    m = 0
    for row_idx in 1:nrows
        has_nonzero = false
        for shelf_idx in 1:p
            qty_col = 3 * shelf_idx
            if dfs[row_idx, qty_col] > 0
                has_nonzero = true
                break
            end
        end
        if has_nonzero
            m += 1
        end
    end
    cost = α * max_cum + β * m
    return max_cum, m, cost
end

"""
    neighbor_ppartition(shelves)::Vector{Vector}

Generate a neighbor solution by relocating one subjob to a different position.

# Arguments
- `shelves`: Current partition as vector of shelves

# Returns
- `Vector{Vector}`: New partition with one subjob relocated

# Algorithm
1. Select a random subjob from any shelf
2. Select a random target shelf (can be same or different)
3. Select a random position in target shelf
4. Remove subjob from source and insert at target position
5. Handle position adjustment if source and target are the same shelf

# Neighborhood Structure
This defines the "move" neighborhood for Simulated Annealing:
- Each move relocates exactly one subjob
- Preserves all job assignments and mold quantities
- Explores different shelf orderings and job distributions

# Edge Cases
- Empty partition returns unmodified copy
- Can move subjob within the same shelf (reordering)

# See Also
- [`neighbor_quantity_redistribution`](@ref) for quantity-based moves
- [`two_step_neighbor`](@ref) for combined neighborhood
"""
function neighbor_ppartition(shelves)
    new_shelves = deepcopy(shelves)
    p = length(new_shelves)
    all_subjobs = []
    for (shelf_idx, shelf) in enumerate(new_shelves)
        for (subjob_idx, subjob) in enumerate(shelf)
            push!(all_subjobs, (shelf_idx=shelf_idx, subjob_idx=subjob_idx, subjob=subjob))
        end
    end
    if isempty(all_subjobs)
        return new_shelves
    end
    selected = rand(all_subjobs)
    source_shelf_idx = selected.shelf_idx
    source_subjob_idx = selected.subjob_idx
    selected_subjob = selected.subjob
    target_shelf_idx = rand(1:p)
    target_shelf = new_shelves[target_shelf_idx]
    max_position = length(target_shelf) + 1
    target_position = rand(1:max_position)
    deleteat!(new_shelves[source_shelf_idx], source_subjob_idx)
    if target_shelf_idx == source_shelf_idx && target_position > source_subjob_idx
        target_position -= 1
    end
    insert!(new_shelves[target_shelf_idx], target_position, selected_subjob)
    return new_shelves
end

"""
    two_step_neighbor(shelves, g, o)::Vector{Vector}

Generate a neighbor solution using adaptive two-step neighborhood strategy.

# Arguments
- `shelves`: Current partition
- `g`: Vector of job identifiers
- `o`: Vector of molds per job (parallel to `g`)

# Returns
- `Vector{Vector}`: New partition created through one or two neighborhood moves

# Algorithm
Adaptive strategy based on selected job type:
1. Select a random subjob
2. Determine job's mold count:
   - **Single-mold job**: Apply two consecutive `neighbor_ppartition` moves
   - **Multi-mold job**: Apply `neighbor_quantity_redistribution` (redistributes quantities across molds)

# Rationale
- Single-mold jobs benefit from repositioning (two moves for better exploration)
- Multi-mold jobs benefit from quantity redistribution across molds
- Provides diverse neighborhood exploration

# Use Cases
Primary neighborhood operator for Simulated Annealing, combining:
- Position-based moves (via `neighbor_ppartition`)
- Quantity-based moves (via `neighbor_quantity_redistribution`)

# See Also
- [`neighbor_ppartition`](@ref) for position-based moves
- [`neighbor_quantity_redistribution`](@ref) for quantity redistribution
"""
function two_step_neighbor(shelves, g, o)
    all_subjobs = []
    for (shelf_idx, shelf) in enumerate(shelves)
        for (subjob_idx, subjob) in enumerate(shelf)
            push!(all_subjobs, subjob)
        end
    end
    if isempty(all_subjobs)
        return deepcopy(shelves)
    end
    selected_subjob = rand(all_subjobs)
    selected_jobid = selected_subjob.jobid
    job_idx = findfirst(==(selected_jobid), g)
    if isnothing(job_idx)
        @warn "Job ID $selected_jobid not found in g vector"
        return deepcopy(shelves)
    end
    num_molds = o[job_idx]
    if num_molds == 1
        intermediate_shelves = neighbor_ppartition(shelves)
    else
        final_shelves = neighbor_quantity_redistribution(shelves, g, zeros(Int, length(g)), o)
        return final_shelves
    end
    final_shelves = neighbor_ppartition(intermediate_shelves)
    return final_shelves
end

"""
    neighbor_quantity_redistribution(shelves, g, n, o)::Vector{Vector}

Generate a neighbor by randomly redistributing quantities across a multi-mold job's molds.

# Arguments
- `shelves`: Current partition
- `g`: Vector of job identifiers
- `n`: Vector of job quantities (currently unused but kept for interface consistency)
- `o`: Vector of molds per job (parallel to `g`)

# Returns
- `Vector{Vector}`: New partition with redistributed mold quantities for one job

# Algorithm
1. Select a random multi-mold job (where o[job] > 1)
2. Locate all subjobs (molds) belonging to this job
3. Calculate total quantity across all molds
4. Generate new random quantity distribution using [`split_quantity_randomly`](@ref)
5. Assign new quantities to molds while preserving their locations

# Neighborhood Structure
This defines the "quantity redistribution" neighborhood:
- Total job quantity remains constant
- Mold positions on shelves remain unchanged
- Only quantity distribution across molds changes

# Invariants
- Total quantity per job preserved
- Number of molds per job preserved
- Shelf assignments of molds preserved

# Edge Cases
- If no multi-mold jobs exist, returns unmodified partition
- Validates that all molds are found before redistribution

# See Also
- [`neighbor_ppartition`](@ref) for position-based moves
- [`two_step_neighbor`](@ref) for combined strategy
"""
function neighbor_quantity_redistribution(shelves, g, n, o)
    new_shelves = deepcopy(shelves)
    multi_mold_jobs = [jobid for (idx, jobid) in enumerate(g) if o[idx] > 1]
    if isempty(multi_mold_jobs)
        return new_shelves
    end
    selected_job = rand(multi_mold_jobs)
    job_idx = findfirst(==(selected_job), g)
    num_molds = o[job_idx]
    mold_locations = []
    for (shelf_idx, shelf) in enumerate(new_shelves)
        for (subjob_idx, subjob) in enumerate(shelf)
            if subjob.jobid == selected_job
                push!(mold_locations, (shelf_idx=shelf_idx, subjob_idx=subjob_idx, 
                                       moldid=subjob.moldid, quantity=subjob.quantity))
            end
        end
    end
    if length(mold_locations) != num_molds
        @warn "Job $selected_job expected $num_molds molds but found $(length(mold_locations))"
        return new_shelves
    end
    total_qty = sum(loc.quantity for loc in mold_locations)
    new_quantities = split_quantity_randomly(total_qty, num_molds)
    sort!(mold_locations, by=x -> x.moldid)
    for (i, loc) in enumerate(mold_locations)
        subjob = new_shelves[loc.shelf_idx][loc.subjob_idx]
        new_shelves[loc.shelf_idx][loc.subjob_idx] = 
            (jobid=subjob.jobid, moldid=subjob.moldid, quantity=new_quantities[i])
    end
    return new_shelves
end

"""
    run(order_dict; log_every::Int=10)::Vector{Vector}

Execute Simulated Annealing algorithm for parallel machine scheduling.

# Arguments
- `order_dict::Dict`: Problem instance dictionary with required keys:
  - `:g` - Job identifiers
  - `:n` - Job quantities
  - `:o` - Molds per job
  - `:p` - Number of parallel shelves/machines
  - `:α` - Shelf load penalty coefficient
  - `:β` - Shelf level penalty coefficient
  - `:Nit` - Maximum iterations
  - `:T0` - Initial temperature
  - `:Tf` - Final (minimum) temperature
  - `:Tj` - Temperature update frequency (iterations between cooling steps)
- `log_every::Int=10`: Logging frequency (iterations per log message, 0 to disable)

# Returns
- `Vector{Vector}`: Best partition found (shelves with assigned subjobs)

# Algorithm: Simulated Annealing
Classical SA with geometric cooling schedule:

1. **Initialization**: Generate random partition
2. **Main Loop**: Until max iterations or min temperature:
   - Generate neighbor using [`two_step_neighbor`](@ref)
   - Calculate cost difference Δ
   - Accept move if:
     - Δ ≤ 0 (improvement), or
     - rand() < exp(-Δ/T) (Metropolis criterion)
   - Update best solution if current improves
   - Cool temperature: T ← T × γ^Tj every Tj iterations
3. **Return**: Best solution found

# Cooling Schedule
Geometric cooling: γ = (Tf / T0)^(1 / Nit)
- Applied every `Tj` iterations
- Ensures temperature reaches `Tf` by iteration `Nit`

# Acceptance Probability
P(accept) = exp(-Δ/T) for worsening moves, allowing escape from local optima

# Output
Logs progress at specified intervals and displays final results including:
- Best cost, shelf levels, and maximum load
- Detailed partition structure
- DataFrames with partition and slot representations

# Examples
```julia
# Standard run with default logging
shelves = run(order_dict)

# Silent mode
shelves = run(order_dict; log_every=0)

# Verbose mode
shelves = run(order_dict; log_every=1)
```

# See Also
- [`two_step_neighbor`](@ref) for neighborhood generation
- [`get_cost_from_shelves`](@ref) for objective evaluation
"""
function run(order_dict; log_every::Int=10)
    @info "Running pure simulation with order_dict:"
    @info order_dict
    @unpack g, n, o, p, α, β, Nit, T0, Tf, Tj = order_dict
    best_shelves = create_random_partition(g, n, o, p)
    current_shelves = deepcopy(best_shelves)
    current_max_sum , current_m, current_cost = get_cost_from_shelves(best_shelves, α, β)
    best_max_sum, best_m, best_cost = current_max_sum, current_m, current_cost
    T = T0
    min_temperature = Tf
    max_iterations = Nit
    γ = (Tf / T0)^(1 / Nit)
    iteration = 0
    accepted = 0
    while iteration < max_iterations && T > min_temperature
        cand_shelves = two_step_neighbor(current_shelves, g, o)
        cand_max_sum , cand_m, cand_cost = get_cost_from_shelves(cand_shelves, α, β)
        Δ = cand_cost - current_cost
        accept = Δ <= 0 || (rand() < exp(-Δ / T))
        if accept
            current_shelves = cand_shelves
            current_max_sum, current_m, current_cost = cand_max_sum, cand_m, cand_cost
            accepted += 1
            if current_cost < best_cost
                best_max_sum, best_m, best_cost = current_max_sum, current_m, current_cost
                best_shelves = deepcopy(cand_shelves)
            end
        end
        if iteration % Tj == 0
            T *= γ^Tj
        end
        iteration += 1
        if log_every > 0 && iteration % log_every == 0
            @info @sprintf(
                "it=%5d T=%.6f cand=%.6f curr=%.6f best=%.6f acc=%4d",
                iteration, T, cand_cost, current_cost, best_cost, accepted)
            accepted = 0
        end
    end
    @info @sprintf("SA finished: it=%d best_cost=%.6f best_m=%d best_max_sum=%.6f", 
                   iteration, best_cost, best_m, best_max_sum)
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

end
