#= Copyright (C) 2024
Nuno David Lopes.
Created:  2024/04/09
Last changed - N. Lopes:2025/03/28 15:54:54
=#

"""
    SimulatedAnnealing

A module that implements a Simulated Annealing optimization algorithm for job scheduling and partitioning problems.

This module includes various functions for building optimization models, solving them, partitioning jobs, and performing the simulated annealing process. It provides functionality for handling different job instances, computing costs, and performing heuristic optimizations across multiple stages. The module supports parallel processing of instances, job transfers between stages, and logging of detailed information for each step of the optimization process.

# Key Functions:
- `model_builder`: Builds the optimization model for a given instance of the job scheduling problem.
- `solver`: Solves the optimization model and computes the objective value.
- `determine_transfer_indexes`: Determines which jobs need to be transferred between stages.
- `solve_and_store`: Solves a job scheduling problem and stores the results.
- `process_instances_group`: Processes multiple job instances, solving them and aggregating the results.
- `transfer_reminder_jobs!`: Transfers jobs to another stage based on certain conditions.
- `move_job!`: Moves a job from one partition to another.
- `simulated_annealing`: Performs the Simulated Annealing optimization process, iterating over job partitions and selecting the best configurations.
- `run_sim`: Executes the full simulation, from initialization to the completion of the simulated annealing process, logging key information along the way.

# Example Usage:
```julia
using SimulatedAnnealing

# Define job instances and settings
instances = ...
order_dict = ...

# Run the simulated annealing optimization process
run_sim(order_file="settings_file.jl")

"""
module SimulatedAnnealing
export run_sim

using DrWatson
@quickactivate "ShoeOptSetupTime"
using Logging, Dates, DataFrames, Gurobi, Random, JuMP, Combinatorics, Printf
include(srcdir("loggers.jl"))



"""
    model_builder(instance; x_cnst_indxs=[], order_dict=order_dict) -> Model

Builds a Gurobi optimization model for job scheduling.

# Arguments
- `instance`: Job parameters (`g`, `n`, `o`).
- `x_cnst_indxs` (optional): Fixed variable constraints.
- `order_dict` (optional): Optimization parameters (`α`, `β`, `p`, `Tl`).

# Returns
- A JuMP model minimizing tardiness and sequencing penalties.

# Example
```julia
model = model_builder(instance; order_dict=order_dict)
optimize!(model)

"""
function model_builder(instance; x_cnst_indxs=[], order_dict=order_dict)
    @debug " Optimization model builder "
    @unpack g, n, o = instance
    @unpack α, β, p, Tl = order_dict

    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "TimeLimit", Tl)
    set_optimizer_attribute(model, "OutputFlag", true)

    set_string_names_on_creation(model, false)

    #set_silent(model)
    unregister(model, :x) # xᵢⱼ
    unregister(model, :t) # tᵢⱼ
    unregister(model, :tbar) #tbarᵢ
    unregister(model, :y) # counter for m

    NJ = length(n)
    NsubJ = sum(o)

    Is = [i for i in 1:NsubJ]
    Js = [j for j in 1:NJ]
    Ks = [[i for i ∈ 1:o[j]] for j ∈ 1:NJ]

    @variable(model, 0 ≤ x[i ∈ Is, j ∈ Js, k ∈ Ks[j]] ≤ 1, Int)
    @variable(model, t[i ∈ Is, j ∈ Js, k ∈ Ks[j]], Int)
    @variable(model, tbar[i ∈ Is])
    @variable(model, y[i ∈ Is])

    @objective(model, Min, sum(α * tbar[i] + β * y[i] for i in Is))
    @constraint(model, [j ∈ Js], sum(t[i, j, k] for i ∈ Is, k ∈ Ks[j]) == n[j])
    @constraint(model, [i ∈ Is, j ∈ Js, k ∈ Ks[j]], tbar[i] ≥ t[i, j, k])
    @constraint(model, [i in Is], sum(x[i, j, k] for j in Js, k ∈ Ks[j]) ≤ p)
    @constraint(model, [i ∈ Is, j ∈ Js, k ∈ Ks[j]], t[i, j, k] ≤ n[j] * x[i, j, k])
    @constraint(model, [i ∈ Is, j ∈ Js, k ∈ Ks[j]], t[i, j, k] ≥ x[i, j, k])
    @constraint(model,
        [i ∈ Is, j ∈ Js, l ∈ Js, u ∈ Ks[j], v ∈ Ks[l]],
        t[i, j, u] ≤ t[i, l, v] + n[j] * (1 - x[i, l, v]))
    @constraint(model,
        [j ∈ Js, k ∈ Ks[j], i ∈ Is[1:end-2], q ∈ Is[i+2:end], l ∈ Is[i+1:q-1]],
        x[l, j, k] ≥ x[i, j, k] + x[q, j, k] - 1)
    @constraint(model, [i ∈ Is, j ∈ Js, k ∈ Ks[j]], y[i] ≥ x[i, j, k])
    @constraint(model, [i ∈ Is[1:end-1]], y[i] ≥ y[i+1])

    @constraint(model, [(j, k) ∈ x_cnst_indxs], x[1, j, k] == 1.0)

    @debug "Model built"

    return model
end


"""
    solver(model, instance; order_dict=order_dict) -> (table_dict, sol_dict)

Solves the optimization model and extracts key results.

# Arguments
- `model`: JuMP model to solve.
- `instance`: Problem data (`n`, `o`, `g`).
- `order_dict` (optional): Optimization parameters (`α`, `β`, `p`).

# Returns
- `table_dict`: Summary metrics (cost, memory, runtime, status).
- `sol_dict`: Solution details (decision variables, indices).

# Example
```julia
table_dict, sol_dict = solver(model, instance)

"""
function solver(model, instance; order_dict=order_dict)
    @debug "Solver for Instance:  $instance "
    @unpack n, o, g = instance
    @unpack α, β, p = order_dict

    optimize!(model)
    MiB = round(Sys.maxrss() / 1048576; digits=1)


    if !is_solved_and_feasible(model; dual=false)
        if !has_values(model)
            error("Model Unsolved: End")
        else
            @warn """
              The model was not solved correctly:
              termination_status : $(termination_status(model))
              primal_status      : $(primal_status(model))
              dual_status        : $(dual_status(model))
              raw_status         : $(raw_status(model))
              """
        end
    end

    stime = round(solve_time(model); digits=10)
    cost = objective_value(model)

    @debug "  objective value = ", cost
    m = sum(value.(model[:y]))
    cost_p = (m > 1) ? α * sum(value.(model[:tbar][1:m-1])) + β * (m - 1) : 0
    @debug "  m value = ", m

    NsubJ = sum(o)
    tf = zeros((NsubJ, NsubJ))
    for i ∈ 1:NsubJ
        aux = 0
        for (j, oj) ∈ enumerate(o)
            col = j + aux
            for k ∈ 1:oj
                aux += (k > 1)
                tf[i, col+k-1] = value(model[:t][i, j, k])
            end
        end
    end
    Dtf = DataFrame(tf, ["Job" * string(g[j]) * "," * string(k) for j ∈ 1:length(n) for k ∈ 1:o[j]])

    xf = zeros((NsubJ, NsubJ))
    for i ∈ 1:NsubJ
        aux = 0
        for (j, oj) ∈ enumerate(o)
            col = j + aux
            for k ∈ 1:oj
                aux += (k > 1)
                xf[i, col+k-1] = value(model[:x][i, j, k])
            end
        end
    end
    Dxf = DataFrame(xf, ["Job" * string(g[j]) * "," * string(k) for j ∈ 1:length(n) for k ∈ 1:o[j]])

    Idxs = [(g[j], k) for j ∈ 1:length(n) for k ∈ 1:o[j]]

    status = termination_status(model)
    table_dict = @strdict α β p m cost cost_p stime MiB status
    sol_dict = @dict Dxf Dtf Idxs m

    @debug "Solver done"

    return table_dict, sol_dict
end


"""
    determine_transfer_indexes(df, sols_dict, i2, step, nshelves) -> Vector

Identifies jobs to transfer based on the last scheduled job.

# Arguments
- `df`: DataFrame containing solution details.
- `sols_dict`: Dictionary of solution variables.
- `i2`: Destination instance.
- `step`: Current step in the process.
- `nshelves`: Shelf capacity constraint.

# Returns
- A list of job indices to transfer, or an empty list if no transfer is needed.

# Example
```julia
transfer_idxs = determine_transfer_indexes(df, sols_dict, i2, step, nshelves)

"""
function determine_transfer_indexes(df, sols_dict, i2, step, nshelves)
    @unpack p = order_dict
    m = Int(df[step, :m])
    last_line_x = collect(values(sols_dict[step][:Dxf][m, :]))
    last_line_t = collect(values(sols_dict[step][:Dtf][m, :]))

    if sum(last_line_x) < nshelves
        val = last_line_t[findfirst(>(0), last_line_t)]
        ixs = findall(>(0), last_line_t)
        indexes = [sols_dict[step][:Idxs][ix] for ix in ixs]

        @debug "Indexes to transfer" indexes
        return transfer_reminder_jobs!(i2, indexes, val)
    end
    return []
end


"""
    solve_and_store(instance, df, sols_dict, x_cnst_indxs=[]; order_dict=order_dict) -> (df, sols_dict)

Builds, solves, and stores the optimization model results.

# Arguments
- `instance`: Problem data for optimization.
- `df`: DataFrame to store solution summaries.
- `sols_dict`: Dictionary to store solution details.
- `x_cnst_indxs` (optional): Indices for constraints.
- `order_dict` (optional): Optimization parameters.

# Returns
- Updated `df` and `sols_dict` with the solution.

# Example
```julia
df, sols_dict = solve_and_store(instance, df, sols_dict)

"""
function solve_and_store(instance, df, sols_dict, x_cnst_indxs=[]; order_dict=order_dict)
    model = model_builder(instance; x_cnst_indxs=x_cnst_indxs, order_dict=order_dict)
    table_dict, sol_dict = solver(model, instance; order_dict=order_dict)

    if table_dict !== nothing
        push!(df, only(DataFrame(table_dict)))
        push!(sols_dict, sol_dict)
    end

    @debug " Instance solved and stored"
    return df, sols_dict
end



"""
    process_instances_group(instances_group; order_dict=order_dict) -> (costs, sols_dict_group)

Solves and processes a group of instance partitions while tracking costs and transfers.

# Arguments
- `instances_group`: A collection of grouped problem instances.
- `order_dict` (optional): Parameters for optimization.

# Returns
- `costs`: List of computed costs for each instance group.
- `sols_dict_group`: List of solution dictionaries for each group.

# Example
```julia
costs, sols_dict_group = process_instances_group(instances_group)

"""
function process_instances_group(instances_group; order_dict=order_dict)
    @unpack p, Pg = order_dict

    costs = []
    sols_dict_group = []
    for i in eachindex(instances_group)
        df = DataFrame(:α => [], :β => [], :p => [], :m => [], :cost => [], :cost_p => [], :stime => [], :MiB => [], :status => [])
        sols_dict = []

        inst = []
        for j ∈ 1:Pg
            push!(inst, deepcopy(instances_group[i][j]))
        end


        # STEP 1
        df, sols_dict = solve_and_store(inst[1], df, sols_dict; order_dict=order_dict)
       
        if Pg == 1
            push!(sols_dict_group, sols_dict)
            push!(costs, df[1, "cost"])
            @debug "Solution Summary $df"
            break
        end

        # STEPS 2:Pg
        transfer_indxs = []
        for j ∈ 2:Pg
            x_indxs = determine_transfer_indexes(df, sols_dict, inst[j], j - 1, p)
            push!(transfer_indxs, x_indxs)
            df, sols_dict = solve_and_store(inst[j], df, sols_dict, x_indxs; order_dict=order_dict)
        end
        push!(sols_dict_group, sols_dict)



        total_cost = 0
        transferred_stages = []

        for j in 1:Pg
            if j == Pg  # Last stage always uses "cost"
                total_cost += df[j, "cost"]
            elseif !isempty(transfer_indxs[j])  # Jobs were transferred to the next stage
                total_cost += df[j, "cost_p"]
                push!(transferred_stages, j + 1)  # Track where jobs were transferred
            else  # No jobs transferred
                total_cost += df[j, "cost"]
            end
        end

        push!(costs, total_cost)

        # Debugging messages
        if !isempty(transferred_stages)
            @debug "Jobs transferred to stages: $(transferred_stages)"
        end
        @debug "Solution Summary $df"

    end

    return costs, sols_dict_group
end


"""
    transfer_reminder_jobs!(dest, indexes, val) -> x_indxs

Transfers remaining jobs to a destination instance and updates its attributes.

# Arguments
- `dest`: Destination instance where jobs are transferred.
- `indexes`: List of job indexes to be transferred.
- `val`: Value associated with the transferred jobs.

# Returns
- `x_indxs`: List of updated indexes after transfer.

# Example
```julia
x_indxs = transfer_reminder_jobs!(dest_instance, job_indexes, transfer_value)

"""
function transfer_reminder_jobs!(dest, indexes, val)
    x_indxs = []
    pos = length(dest[:n])

    g_ids = [idx[1] for idx ∈ indexes]
    count_g_ids = [count(==(elmt), g_ids) for elmt ∈ unique(g_ids)]
    unique!(g_ids)

    values = [val * count_g_ids[i] for i ∈ eachindex(g_ids)]

    for (g_id, count_g_id, value) ∈ zip(g_ids, count_g_ids, values)
        push!(dest[:g], g_id)
        push!(dest[:o], count_g_id)
        push!(dest[:n], value)
    end

    for i ∈ eachindex(g_ids)
        append!(x_indxs, [(pos + i, j) for j in 1:dest[:o][pos+i]])
    end

    @debug "Indexes transfered" x_indxs
    return x_indxs
end



"""
    move_job!(src, dest, jobindex)

Moves a job from the source instance (`src`) to the destination instance (`dest`) at the specified job index.

# Arguments
- `src`: The source instance from which the job will be moved.
- `dest`: The destination instance to which the job will be moved.
- `jobindex`: The index of the job to be moved.

# Example
```julia
move_job!(source_instance, dest_instance, job_index)

"""
function move_job!(src, dest, jobindex)
    for key in [:g, :n, :o]
        push!(dest[key], src[key][jobindex])
        deleteat!(src[key], jobindex)
    end
    @debug "Job moved"
end



"""
    partition_jobs(n, o, g, p, Pg)

Partitions jobs into `Pg` groups such that the total size of each group is at least `p`.

# Arguments
- `n`: A vector of job identifiers.
- `o`: A vector of job sizes corresponding to each job in `n`.
- `g`: A vector of job groups corresponding to each job in `n`.
- `p`: The minimum size each group should have.
- `Pg`: The number of groups to partition the jobs into.

# Returns
- `ns`: A list of job identifiers partitioned into `Pg` groups.
- `os`: A list of job sizes partitioned into `Pg` groups.
- `gs`: A list of job groups partitioned into `Pg` groups.

# Example
```julia
ns, os, gs = partition_jobs(job_ids, job_sizes, job_groups, min_size, num_groups)

"""
function partition_jobs(n, o, g, p, Pg)
    while true
        ns, os, gs = [Vector() for _ ∈ 1:Pg], [Vector() for _ ∈ 1:Pg], [Vector() for _ ∈ 1:Pg]
        sizes = zeros(Pg)

        for j in eachindex(n)
            idx = rand(1:Pg)
            push!(ns[idx], n[j])
            push!(os[idx], o[j])
            push!(gs[idx], g[j])
            sizes[idx] += o[j]
        end

        if all(s -> s ≥ p, sizes)
            @debug "Partitioned jobs in $Pg groups"
            return ns, os, gs
        end
    end
end


"""
    generate_instance_groups(instances, Pg)

Generates all permutations of the instances based on the given partition size `Pg`. The function only supports partition sizes between 1 and 3.

# Arguments
- `instances`: A vector of instances to be partitioned.
- `Pg`: The number of groups to partition the instances into (must be between 1 and 3).

# Returns
- A list of permutations of the instances for the given partition size `Pg`. 
  If `Pg` is not between 1 and 3, an error is raised.

# Example
```julia
groups = generate_instance_groups([inst1, inst2, inst3], 2)
"""
function generate_instance_groups(instances, Pg)
    if 1 ≤ Pg ≤ 3
        return collect(permutations(instances[1:Pg]))
    else
        error("Partition size not supported: EXIT")
        return []
    end
end


"""
    simulated_annealing_init(; order_dict=order_dict)

Initializes the simulated annealing process by partitioning jobs, generating instance groups, and processing them to find the best solution.

# Arguments
- `order_dict`: A dictionary containing various parameters needed for the initialization, such as:
  - `g`: A vector of group sizes.
  - `n`: A vector of job sizes.
  - `o`: A vector of job orders.
  - `p`: A threshold parameter.
  - `α`, `β`: Cost coefficients.
  - `Nit`: Number of iterations (though not used in this function).
  - `Pg`: Number of groups for partitioning jobs.

# Returns
- `instances`: The selected best instance group based on the simulated annealing initialization process.
- `cost`: The cost associated with the best instance group.
- `sols_dict`: The solution dictionary corresponding to the best instance group.

# Description
This function partitions jobs into `Pg` groups and generates permutations of instances. It processes these instances to compute costs, and selects the one with the minimum cost as the best initial solution for the simulated annealing algorithm.

# Example
```julia
instances, cost, sols_dict = simulated_annealing_init(order_dict=some_dict)

"""
function simulated_annealing_init(; order_dict=order_dict)
    @unpack g, n, o, p, α, β, Nit, Pg = order_dict
    @debug "Simulated Annealing Initialization Function"
    local instances = []

    ns, os, gs = partition_jobs(n, o, g, p, Pg)


    for (n, o, g) ∈ zip(ns, os, gs)
        instance = @dict α β n o g p
        push!(instances, instance)
    end

    instances_group = generate_instance_groups(instances, Pg)
    @info "Instances group length $(length(instances_group))"
    
    costs, sols_dict_group = process_instances_group(instances_group)

    #@info "Sols_dict_group $sols_dict_group" 
    @info "Costs  $costs"
    cost, min_indx = findmin(costs)
    @info "Best cost $cost at instances $min_indx"
    instances = instances_group[min_indx]
    @debug "Instances  $instances"
    sols_dict = sols_dict_group[min_indx]

    for (i,sols) ∈ enumerate(sols_dict_group)
        @debug "Instance permutation $i"
        @debug sols
    end

    return instances, cost, sols_dict

end


"""
    neighbour_partition!(instances, g, p, Pg)

Performs a job transfer between partitions to improve the current solution, selecting a random job and attempting to move it to another partition if certain conditions are met.

# Arguments
- `instances`: A list of job partitions.
- `g`: A vector of job IDs to choose from.
- `p`: A threshold parameter used to check if the job can be moved.
- `Pg`: The number of partitions (must be at least 2).

# Returns
- The `jobid` of the job that was successfully moved, or raises an error if no suitable job is found after a maximum number of attempts.

# Description
This function attempts to move a randomly chosen job from one partition to another, ensuring that the destination partition meets certain constraints (e.g., the sum of orders in the source partition must be greater than or equal to `p`). It tries a maximum of `max_counter` iterations, and will exit with an error if no valid job transfer is found.

# Example
```julia
jobid = neighbour_partition!(instances, g, p, Pg)

"""
function neighbour_partition!(instances, g, p, Pg)

    if Pg < 2 || Pg > length(instances)
        error("Partition size not supported: EXIT")
    end

    counter = 0
    max_counter = length(g)

    while counter < max_counter
        counter += 1
        @debug "Counter $counter for neighbour_partition"
        jobid = rand(g)
        for partition ∈ collect(permutations(instances[1:Pg]))
            src = partition[1]
            dest = rand(partition[2:end])  # Randomly pick one of the remaining destinations

            if jobid ∈ src[:g] && sum(src[:o]) ≥ p
                move_job!(src, dest, findfirst(==(jobid), src[:g]))
                @debug "neighbour_partition done: job $jobid moved"
                return jobid
            end
        end
        if counter == max_counter
            error("No available neighbours: EXIT")
        end
    end
end


"""
    simulated_annealing(input_instances, input_cost, input_sols_dict; order_dict=order_dict)

Performs the Simulated Annealing optimization algorithm to minimize the cost of a set of instances over a number of iterations.

# Arguments
- `input_instances`: The initial set of job instances.
- `input_cost`: The initial cost of the current solution.
- `input_sols_dict`: The initial solution dictionary.
- `order_dict`: A dictionary containing various parameters:
  - `g`: Group identifiers.
  - `p`: Partition size.
  - `Nit`: Number of iterations.
  - `Pg`: Number of partitions.
  - `Gl`: Global time limit.
  - `T0`: Initial temperature.
  - `Tf`: Final temperature.
  - `Tj`: Temperature adjustment interval.

# Returns
- `best_instances`: The best solution found after all iterations.
- `best_cost`: The cost of the best solution.
- `best_sols_dict`: The solution dictionary for the best solution.
- `cur_costs_list`: A list of costs recorded at each iteration.

# Description
This function runs the Simulated Annealing algorithm, starting with an initial solution, and iteratively improves the solution by performing random job moves between partitions. The temperature is gradually decreased according to a cooling factor, allowing the algorithm to escape local minima.

The algorithm accepts random job moves, but only accepts worse solutions based on a probabilistic function of the temperature. The temperature decreases as the algorithm progresses, reducing the probability of accepting worse solutions over time.

# Example
```julia
best_instances, best_cost, best_sols_dict, cur_costs_list = simulated_annealing(input_instances, input_cost, input_sols_dict; order_dict=order_dict)

"""
function simulated_annealing(input_instances, input_cost, input_sols_dict; order_dict=order_dict)
    @unpack g, p, Nit, Pg, Gl = order_dict
    @unpack T0, Tf, Tj = order_dict
    @info "Simulated Annealing with $Nit iterations"

    # Cooling factor
    γ = (Tf / T0)^(1 / Nit)

    # Initialize best and current states
    best_instances, best_cost, best_sols_dict = deepcopy(input_instances), input_cost, deepcopy(input_sols_dict)
    cur_instances, cur_cost, cur_sols_dict = deepcopy(best_instances), best_cost, deepcopy(best_sols_dict)
    cur_costs_list = [best_cost]

    Temp = T0
    time_counter = 0

    for i ∈ 2:Nit
        iter_time_start = time()
        jobid = neighbour_partition!(cur_instances, g, p, Pg)
        @debug "Selected job to switch: $jobid"

        instances_group = generate_instance_groups(cur_instances, Pg)
        alt_costs, alt_sols_dict_group = process_instances_group(instances_group)

        @debug "Alternative costs  $alt_costs"

        alt_cost, best_alt_idx = findmin(alt_costs)
        @info "Best cost $cost at instances $best_alt_idx"
        alt_instances, alt_sols_dict = instances_group[best_alt_idx], alt_sols_dict_group[best_alt_idx]

        if i % Tj == 0
            Temp *= γ^Tj # Decrease temperature
        end

        Δ = (cur_cost - alt_cost) / Temp

        formated_cur_cost = @sprintf("%04.1f", cur_cost)
        formated_alt_cost = @sprintf("%04.1f", alt_cost)
        formated_best_cost = @sprintf("%04.1f", best_cost)
        formated_Δ = @sprintf("%+.5e", Δ)
        formated_exp_Δ = @sprintf("%.5e", exp(Δ))
        formated_Temp = @sprintf("%09.6f", Temp)
        formated_Time = @sprintf("%07.1f", time_counter)
        formated_iter = @sprintf("%04i", i)
        @info "Iter $(formated_iter) | cur  $(formated_cur_cost) | alt $(formated_alt_cost) | best $(formated_best_cost) | Temp  $(formated_Temp) | Δ $(formated_Δ) | exp(Δ)  $(formated_exp_Δ)| sa_timer $(formated_Time) |"

        if alt_cost ≤ cur_cost
            cur_cost, cur_instances, cur_sols_dict = alt_cost, deepcopy(alt_instances), deepcopy(alt_sols_dict)
            if alt_cost < best_cost
                best_cost, best_instances, best_sols_dict = alt_cost, deepcopy(alt_instances), deepcopy(alt_sols_dict)
            end
        else

            if exp(Δ) > rand()
                cur_instances, cur_cost, cur_sols_dict = deepcopy(alt_instances), alt_cost, deepcopy(alt_sols_dict)
            end
        end

        push!(cur_costs_list, cur_cost)
        @debug "Iteration $i/$Nit - Best: $best_cost, Current: $cur_cost, Alternative: $alt_cost"

        iter_time_stop = time()
        time_counter += iter_time_stop - iter_time_start
        @debug "Simulated Annealing Elapsed time $time_counter seconds"

        if time_counter > Gl
            @debug "Global Time Exceeded:  $time_counter  >  $Gl : loop break."
            break
        end
    end

    return best_instances, best_cost, best_sols_dict, cur_costs_list
end


"""
    run_sim(; order_file=nothing)

Runs the Simulated Annealing optimization process using a provided order file, logging the process and results.

# Arguments
- `order_file`: The path to a file containing the order settings. If not provided, the function will throw an error.

# Returns
- `nothing`: The function does not return any values but logs the results to a specified location.

# Description
This function initializes the simulation, performs a Simulated Annealing optimization, and logs key information at each step. It begins by loading the settings from the provided order file. If the partition size (`Pg`) is 1, the function skips the Simulated Annealing process.

The function tracks the time spent on initialization, Simulated Annealing, and logs details on the best solutions found, including the `m` values and the heuristic cost.

If the order file is not provided, the function raises an error and exits.

# Example
```julia
run_sim(order_file="path_to_order_file")

"""
function run_sim(; order_file=nothing)
    if order_file == nothing
        error("No order file provided: EXIT")
        return nothing
    end

    include(datadir("settings", order_file))
    prefix = datadir("sims", order_file[1:end-3] * order_dict[:Oid])

    with_logger(demux_logger(prefix)) do
        elapsed_time = @elapsed begin
            @info "Simulation date  $(now())"
            @info "Order Dictionary $(order_dict)"

            init_instances, init_cost, init_sols_dict = simulated_annealing_init(; order_dict=order_dict)

            @info "Simulated Annealing Init"
            @info "Initial instances: $init_instances"
            @info "Initial cost: $init_cost"
            @info "Initial Solutions (ts) and m values:"
            for sol ∈ init_sols_dict
                @info sol[:Dtf]
                @info "m value = $(sol[:m])"
            end
        end
        
        @unpack Pg = order_dict
        if Pg == 1
            @info "Partition size is 1, no need to run simulated annealing"
            @info "Elapsed time $(elapsed_time) seconds"
            return nothing
        end
        @info "Elapsed time $(elapsed_time) seconds"

        elapsed_time += @elapsed begin
            best_instances, best_cost, best_sols_dict, cur_costs_list = simulated_annealing(init_instances, init_cost, init_sols_dict; order_dict=order_dict)
            @info "Simulated Annealing End"

            # Info on best solutions
            @info "Best Instances: $best_instances"
            @info "Best Solutions (ts):"

            h_m = 0
            @unpack p = order_dict
            for (i, sol) ∈ enumerate(best_sols_dict)
                @info sol[:Dtf]
                m = Int(sol[:m])
                last_line_x = collect(values(sol[:Dxf][m, :]))
                if sum(last_line_x) < p && i < length(best_sols_dict)
                    h_m += sol[:m] - 1
                else
                    h_m += sol[:m]
                end
            end
            @info "Heuristic solutions Best cost: $best_cost"
            @info "Heuristic solution  m value: $(h_m)"

            @info "Heuristic iterations Current costs: $cur_costs_list"
        end
        @info "Elapsed time $(elapsed_time) seconds"
    end
end

# End of module
end

