#=  Copyright (C) 2024
Nuno David Lopes.
Created:  2024/04/09
Last changed - N. Lopes: 2025/11/26 09:47:15
=#

module SimulatedAnnealing
export run, create_random_partition, partition_to_dataframe
export extract_slots_from_partition, neighbor_ppartition, neighbor_quantity_redistribution
export get_cost_from_shelves, get_cost_from_slots, split_quantity_randomly, two_step_neighbor

using DrWatson
@quickactivate "ShoeOptSetupTime"

using Logging, Dates, DataFrames, Random, Combinatorics, Printf


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
