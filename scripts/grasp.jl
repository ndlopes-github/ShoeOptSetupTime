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
    even_split(total::Int, parts::Int) :: Vector{Int}

Deterministically split `total` into `parts` non-negative integers that sum to `total`.
The values are as equal as possible; earlier parts get the extra 1 when `total % parts != 0`.
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
    create_all_subjobs(g, n, o; split_method = :random)

Create the list of subjobs (jobid, moldid, quantity) for all jobs.
split_method:
- :random  -> uses SimulatedAnnealing.split_quantity_randomly (current behavior)
- :even    -> uses deterministic even_split to remove randomness
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
    select_subjob(subjobs; mode=:probabilistic)

Select a subjob from the list.
mode:
- :probabilistic -> roulette-wheel proportional to quantity (current behavior)
- :greedy       -> deterministically pick the subjob with max quantity; ties resolved by first occurrence
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

function get_shelf_loads(shelves)
    loads = zeros(Int, length(shelves))
    for (shelf_idx, shelf) in enumerate(shelves)
        if !isempty(shelf)
            loads[shelf_idx] = sum(sj.quantity for sj in shelf)
        end
    end
    return loads
end

function place_subjob_greedy(subjob, shelves)
    new_shelves = deepcopy(shelves)
    loads = get_shelf_loads(new_shelves)
    new_loads = loads .+ subjob.quantity
    min_load_shelf = argmin(new_loads)
    push!(new_shelves[min_load_shelf], subjob)
    return new_shelves
end

"""
    create_grasp_partition(g, n, o, p; selection=:probabilistic, split_method=:random)

Build one partition using GRASP-style construction.
- selection=:probabilistic uses roulette selection; :greedy is deterministic max-quantity selection.
- split_method=:random uses random mold splits; :even is deterministic.
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
    run(order_dict; log_every=10, selection=:probabilistic, split_method=:random, iterations=nothing)

Run GRASP for `Nit` iterations from `order_dict` (or `iterations` if provided), returning the best shelves.

To obtain a pure greedy solution: set `selection=:greedy`, `split_method=:even` and run 1 iteration
either by passing `iterations=1` or ensuring `order_dict[:Nit] == 1`.
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
    run_greedy(order_dict; log_every=0)

Compute a single deterministic greedy partition (no randomness) and return summary.
This uses selection=:greedy and split_method=:even with one construction.

Returns a NamedTuple with fields: shelves, max_sum, m, cost, time.
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
