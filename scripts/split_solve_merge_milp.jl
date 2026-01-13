#=  Copyright (C) 2024
Nuno David Lopes.
Created:  2024/04/09
Last changed - N. Lopes:2025/12/07 09:55:04
=#

module SplitSolveMergeMILP
export run

using DrWatson
@quickactivate "ShoeOptSetupTime"

using Logging, Dates, DataFrames, Gurobi, HiGHS, Random, JuMP, Combinatorics, Printf

include(srcdir("loggers.jl"))

"""
    model_builder(instance; x_cnst_indxs=[], order_dict)::Model

Build a Mixed-Integer Linear Programming (MILP) model for the parallel machine scheduling problem.

# Arguments
- `instance::Dict`: Instance data containing:
  - `:g` - Job identifiers
  - `:n` - Job quantities
  - `:o` - Number of molds per job
- `x_cnst_indxs::Vector=[]`: Constraint indices to fix specific x variables (used for sequential solving)
- `order_dict::Dict`: Configuration dictionary with:
  - `:α` - Penalty coefficient for maximum shelf load
  - `:β` - Penalty coefficient for number of shelf levels
  - `:p` - Number of parallel shelves (machines)
  - `:Tl` - Time limit for solver (seconds)
  - `:solver_name` - Solver to use ("Gurobi" or "HiGHS")

# Returns
- `Model`: Configured JuMP optimization model ready to be solved

# Model Formulation

## Decision Variables
- `x[i,j,k]` ∈ {0,1}: Binary indicator if subjob (j,k) is assigned to shelf i
- `t[i,j,k]` ∈ ℤ₊: Quantity of job j, mold k assigned to shelf i
- `tbar[i]` ∈ ℝ₊: Maximum cumulative quantity on shelf i
- `y[i]` ∈ {0,1}: Indicator if shelf i is used

## Objective
Minimize: Σᵢ (α × tbar[i] + β × y[i])

## Constraints
1. Job completion: Sum of quantities equals total requirement
2. Shelf capacity: At most p subjobs per shelf
3. Assignment consistency: Quantity constraints based on x
4. Cumulative tracking: tbar captures maximum load
5. Continuity: No gaps in shelf usage
6. Fixed assignments: Optional x_cnst_indxs constraints

# Solver Configuration

## Gurobi
- Time limit, threads, memory management
- Presolve and sparsification enabled

## HiGHS
- Time limit and console logging

# See Also
- [`solver`](@ref) for solving the model
- [`solve_and_store`](@ref) for combined build and solve
"""
function model_builder(instance; x_cnst_indxs=[], order_dict)
    @debug " Optimization model builder "
    @unpack g, n, o = instance
    @unpack α, β, p, Tl, solver_name = order_dict

    # Select optimizer based on solver_name parameter
    if solver_name == "Gurobi"
        model = Model(Gurobi.Optimizer)
        set_optimizer_attribute(model, "TimeLimit", Tl)
        set_optimizer_attribute(model, "OutputFlag", true)
        set_optimizer_attribute(model, "Threads", 30)
        set_optimizer_attribute(model, "NodefileStart", 16.0)
        set_optimizer_attribute(model, "NodefileDir", "/tmp")
        set_optimizer_attribute(model, "Presolve", 1)
        set_optimizer_attribute(model,"SoftMemLimit", 56) #In GB
        set_optimizer_attribute(model, "PreSparsify", 1)
    elseif solver_name == "HiGHS"
        model = Model(HiGHS.Optimizer)
        set_optimizer_attribute(model, "time_limit", Float64(Tl))
        set_optimizer_attribute(model, "log_to_console", true)
    else
        error("Unsupported solver: $solver_name. Use 'Gurobi' or 'HiGHS'")
    end

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
    @variable(model, t[i ∈ Is, j ∈ Js, k ∈ Ks[j]] ≤ sum(n), Int)
    @variable(model, tbar[i ∈ Is])
    @variable(model, y[i ∈ Is])

    @objective(model, Min, sum(α * tbar[i] + β * y[i] for i in Is))
    @constraint(model, [j ∈ Js], sum(t[i, j, k] for i ∈ Is, k ∈ Ks[j]) == n[j])
    @constraint(model, [i ∈ Is, j ∈ Js, k ∈ Ks[j]], tbar[i] ≥ t[i, j, k])
    @constraint(model, [i in Is], sum(x[i, j, k] for j in Js, k ∈ Ks[j]) ≤ p)
    @constraint(model, [i ∈ Is, j ∈ Js, k ∈ Ks[j]], t[i, j, k] ≤ n[j] * x[i, j, k])
    @constraint(model, [i ∈ Is, j ∈ Js, k ∈ Ks[j]], t[i, j, k] ≥ x[i, j, k])


    @constraint(model,
        [i ∈ Is, l ∈ Js, v ∈ Ks[l]],
        tbar[i] ≤ t[i, l, v] + maximum(n) * (1 - x[i, l, v]))

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
	solver(model, instance; order_dict)

Solves the JuMP optimization model and extracts solution details.

This function optimizes the MILP formulation using the configured solver (Gurobi or HiGHS),
extracts decision variables, computes objective values, and formats results for analysis.
It handles both optimal and feasible solutions, providing appropriate warnings for suboptimal terminations.

# Arguments
- `model::JuMP.Model`: The built optimization model with variables and constraints.
- `instance::Dict`: Instance data containing keys `:n` (job quantities), `:o` (molds per job), `:g` (job IDs).
- `order_dict::Dict`: Configuration with keys `:α` (setup cost weight), `:β` (machine cost weight), `:p` (available shelves).

# Returns
- `table_dict::Dict`: Solution summary with keys:
  - `:cost` - Total objective value
  - `:cost_p` - Partial cost excluding last machine's setup
  - `:m` - Number of machines used
  - `:stime` - Solve time in seconds
  - `:MiB` - Memory usage in MiB
  - `:status` - Solver termination status
- `sol_dict::Dict`: Detailed solution with keys:
  - `:Dxf` - DataFrame of assignment variables (machines × job-molds)
  - `:Dtf` - DataFrame of quantity variables (machines × job-molds)
  - `:Idxs` - Tuple array of (job_id, mold) pairs
  - `:m` - Number of machines as integer

# Notes
- Generates warnings if model terminates with suboptimal status
- Memory usage measured via `Sys.maxrss()`
- Uses `@debug` macros for detailed logging during solve process
- See also: [`model_builder`](@ref)

# Examples
```julia
instance = Dict(:n => [100, 80], :o => [2, 3], :g => [1, 2], ...)
model = model_builder(instance; order_dict=order_dict)
table_dict, sol_dict = solver(model, instance; order_dict=order_dict)
println("Machines used: ", sol_dict[:m])
```
"""
function solver(model, instance; order_dict)
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
	determine_transfer_indexes(df, sols_dict, i2, step, nshelves; order_dict)

Determines which jobs should be transferred from the current stage to the next.

In the Split-Solve-Merge algorithm, when a machine in the current stage solution has unfilled
shelves (i.e., fewer than `p` job-mold assignments), those jobs are transferred to the next stage
to continue processing. This function identifies such jobs based on the last machine's assignments.

# Arguments
- `df::DataFrame`: Solution summary table with cost and machine count per stage.
- `sols_dict::Vector`: Detailed solutions for each stage, containing assignment and quantity matrices.
- `i2::Dict`: Next stage instance to which jobs will be transferred.
- `step::Int`: Current stage number (1-indexed).
- `nshelves::Int`: Number of available shelves (`:p` parameter).
- `order_dict::Dict`: Configuration parameters.

# Returns
- `indexes::Vector{Tuple{Int,Int}}`: List of (job_id, mold_number) pairs to transfer, or empty vector if no transfer needed.

# Algorithm
1. Extract last machine (row `m`) from assignment matrix `xᵢⱼ`
2. If `sum(x[m,:]) < p`, collect jobs assigned to that machine
3. Call `transfer_reminder_jobs!` to modify `i2` instance in-place
4. Return constraint indexes for the next stage's model

# Notes
- Transfer only occurs when last machine is underutilized
- See also: [`transfer_reminder_jobs!`](@ref), [`process_instances_group`](@ref)

# Examples
```julia
indexes = determine_transfer_indexes(df, sols_dict, next_instance, 1, 15; order_dict=od)
# Returns: [(1,1), (1,2), (3,1)] if these jobs need transfer
```
"""
function determine_transfer_indexes(df, sols_dict, i2, step, nshelves; order_dict)
    @unpack p = order_dict
    m = round(Int, df[step, :m])
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
	solve_and_store(instance, df, sols_dict, x_cnst_indxs = []; order_dict)

Solves a single instance and stores results in the provided data structures.

This is a convenience wrapper that builds the model, solves it, and appends results to the 
tracking DataFrame and solutions dictionary. It also handles memory cleanup for large-scale runs.

# Arguments
- `instance::Dict`: Instance data with keys `:n`, `:o`, `:g`, `:α`, `:β`, `:p`.
- `df::DataFrame`: Accumulator for solution summaries (modified in-place).
- `sols_dict::Vector`: Accumulator for detailed solutions (modified in-place).
- `x_cnst_indxs::Vector{Tuple{Int,Int}}`: Optional. Assignment constraints from previous stage transfers (default: `[]`).
- `order_dict::Dict`: Configuration parameters.

# Returns
- `df::DataFrame`: Updated with new row containing cost, machines, solve time, status.
- `sols_dict::Vector`: Updated with new solution dictionary containing `Dxf`, `Dtf`, `Idxs`, `m`.

# Side Effects
- Calls [`model_builder`](@ref) to create JuMP model
- Calls [`solver`](@ref) to optimize and extract results
- Explicitly frees model memory via `empty!(model)` for memory efficiency

# Notes
- The `x_cnst_indxs` parameter enforces specific assignments (e.g., transferred jobs must be assigned to machine 1)
- Used extensively in [`process_instances_group`](@ref) for multi-stage solves
- See also: [`model_builder`](@ref), [`solver`](@ref)

# Examples
```julia
df = DataFrame(:cost => [], :m => [], ...)
sols = []
df, sols = solve_and_store(instance, df, sols; order_dict=od)
println("Total cost: ", df[end, :cost])
```
"""
function solve_and_store(instance, df, sols_dict, x_cnst_indxs=[]; order_dict)
    model = model_builder(instance; x_cnst_indxs=x_cnst_indxs, order_dict=order_dict)
    table_dict, sol_dict = solver(model, instance; order_dict=order_dict)

    if table_dict !== nothing
        push!(df, only(DataFrame(table_dict)))
        push!(sols_dict, sol_dict)
    end

    # Free model memory explicitly
    try
        empty!(model)
    catch
        # Older JuMP versions may not support empty!
    end
    model = nothing

    @debug " Instance solved and stored"
    return df, sols_dict
end



"""
	process_instances_group(instances_group; order_dict)

Processes a group of instances representing different stage permutations in the Split-Solve-Merge algorithm.

For each permutation of instances (representing different processing stage orders), this function:
1. Solves the first stage normally
2. For subsequent stages (if `Pg > 1`), determines job transfers and solves with transfer constraints
3. Computes total cost accounting for partial costs when jobs are transferred between stages

The Split-Solve-Merge algorithm partitions jobs into stages and solves each stage sequentially,
transferring underutilized jobs from one stage to the next for optimal resource usage.

# Arguments
- `instances_group::Vector{Vector{Dict}}`: Collection of instance permutations, where each permutation is a sequence of stage instances.
- `order_dict::Dict`: Configuration with keys:
  - `:p` - Available shelves per machine
  - `:Pg` - Number of processing stages (partition groups)

# Returns
- `costs::Vector{Float64}`: Total cost for each permutation in `instances_group`.
- `sols_dict_group::Vector{Vector{Dict}}`: Detailed solutions for each permutation (nested: permutation → stage → solution).

# Algorithm
- **Stage 1**: Solve normally without constraints
- **Stages 2 to Pg**: 
  1. Determine transferred jobs from previous stage
  2. Add assignment constraints for transferred jobs
  3. Solve current stage
- **Cost Calculation**: Use `:cost_p` (partial) if jobs transferred to next stage, otherwise use `:cost`

# Notes
- If `Pg == 1`, only solves a single stage (exact MILP, no merging)
- Transfer decisions based on shelf utilization in last machine of each stage
- See also: [`determine_transfer_indexes`](@ref), [`solve_and_store`](@ref), [`simulated_annealing_init`](@ref)

# Examples
```julia
instances_group = [
    [inst1_perm1, inst2_perm1, inst3_perm1],
    [inst2_perm2, inst1_perm2, inst3_perm2]
]
costs, sols = process_instances_group(instances_group; order_dict=od)
best_cost, best_idx = findmin(costs)
```
"""
function process_instances_group(instances_group; order_dict)
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
            x_indxs = determine_transfer_indexes(df, sols_dict, inst[j], j - 1, p; order_dict=order_dict)
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
	transfer_reminder_jobs!(dest, indexes, val)

Transfers reminder jobs to the destination instance, modifying it in-place.

When a stage's last machine has underutilized shelves, jobs assigned to it are transferred
to the next stage for continued processing. This function adds those jobs to the destination
instance's job lists and returns constraint indexes to enforce their assignment in the next solve.

# Arguments
- `dest::Dict`: Destination instance to be modified. Must have keys `:g`, `:o`, `:n` (all vectors).
- `indexes::Vector{Tuple{Int,Int}}`: List of (job_id, mold_number) pairs to transfer.
- `val::Float64`: Quantity value to assign to each transferred job in `:n`.

# Returns
- `x_indxs::Vector{Tuple{Int,Int}}`: Constraint indexes `(j, k)` for the destination model, ensuring transferred jobs are assigned to machine 1.

# Side Effects
- Modifies `dest[:g]` by appending unique job IDs from `indexes`
- Modifies `dest[:o]` by appending mold counts per transferred job
- Modifies `dest[:n]` by appending quantities (based on `val` and mold counts)

# Algorithm
1. Extract unique job IDs from `indexes`
2. Count mold occurrences per job ID
3. Compute transfer quantities: `val × mold_count`
4. Append to `dest[:g]`, `dest[:o]`, `dest[:n]`
5. Generate constraint indexes for all transferred job-molds

# Notes
- The `!` suffix indicates in-place modification
- Transfer constraints force transferred jobs onto machine 1 in the next stage
- See also: [`determine_transfer_indexes`](@ref), [`process_instances_group`](@ref)

# Examples
```julia
dest_instance = Dict(:g => [1], :o => [2], :n => [50.0])
indexes = [(2,1), (2,2), (3,1)]  # Jobs 2 and 3 with molds
x_indxs = transfer_reminder_jobs!(dest_instance, indexes, 30.0)
# dest_instance now has jobs [1, 2, 3] with updated :o and :n
# x_indxs: [(2,1), (2,2), (3,1)] for constraint enforcement
```
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

Moves a job from source instance to destination instance (in-place).

Transfers all job attributes (`:g`, `:n`, `:o`) from the source instance to the destination,
removing it from the source. Used during neighborhood generation in simulated annealing to
rebalance job allocations across stages.

# Arguments
- `src::Dict`: Source instance with keys `:g`, `:o`, `:n` (all vectors).
- `dest::Dict`: Destination instance with keys `:g`, `:o`, `:n` (all vectors).
- `jobindex::Int`: Index of the job in `src` vectors to move (1-based).

# Side Effects
- Appends `src[:g][jobindex]`, `src[:o][jobindex]`, `src[:n][jobindex]` to `dest` vectors
- Removes elements at `jobindex` from all `src` vectors

# Notes
- The `!` suffix indicates in-place modification
- Used by [`neighbour_partition!`](@ref) during simulated annealing
- Ensures job counts remain consistent: `length(src[:g]) == length(src[:o]) == length(src[:n])`
- See also: [`switch_jobs!`](@ref), [`neighbour_partition!`](@ref)

# Examples
```julia
src = Dict(:g => [1, 2, 3], :o => [2, 3, 1], :n => [100, 80, 50])
dest = Dict(:g => [4], :o => [2], :n => [60])
move_job!(src, dest, 2)  # Move job 2 from src to dest
# src: :g => [1, 3], :o => [2, 1], :n => [100, 50]
# dest: :g => [4, 2], :o => [2, 3], :n => [60, 80]
```
"""
function move_job!(src, dest, jobindex)
    for key in [:g, :n, :o]
        push!(dest[key], src[key][jobindex])
        deleteat!(src[key], jobindex)
    end
    @debug "Job moved"
end



"""
	switch_jobs!(inst1, inst2, jobindex1, jobindex2)

Swaps jobs between two instances (in-place).

Exchanges all job attributes (`:g`, `:n`, `:o`) at specified indices between two instances.
Used during neighborhood generation when direct job moves would violate the minimum shelf
constraint (`sum(o) ≥ p`).

# Arguments
- `inst1::Dict`: First instance with keys `:g`, `:o`, `:n` (all vectors).
- `inst2::Dict`: Second instance with keys `:g`, `:o`, `:n` (all vectors).
- `jobindex1::Int`: Index of the job in `inst1` to swap (1-based).
- `jobindex2::Int`: Index of the job in `inst2` to swap (1-based).

# Side Effects
- Exchanges `inst1[:g][jobindex1] ↔ inst2[:g][jobindex2]`
- Exchanges `inst1[:o][jobindex1] ↔ inst2[:o][jobindex2]`
- Exchanges `inst1[:n][jobindex1] ↔ inst2[:n][jobindex2]`

# Notes
- The `!` suffix indicates in-place modification
- Used by [`neighbour_partition!`](@ref) when `sum(src[:o]) == p` (exact shelf match)
- Swaps preserve job counts in both instances
- Typically swaps jobs with matching mold counts to maintain shelf feasibility
- See also: [`move_job!`](@ref), [`neighbour_partition!`](@ref)

# Examples
```julia
inst1 = Dict(:g => [1, 2], :o => [2, 3], :n => [100, 80])
inst2 = Dict(:g => [3, 4], :o => [2, 2], :n => [60, 50])
switch_jobs!(inst1, inst2, 2, 1)  # Swap job 2 in inst1 with job 3 in inst2
# inst1: :g => [1, 3], :o => [2, 2], :n => [100, 60]
# inst2: :g => [2, 4], :o => [3, 2], :n => [80, 50]
```
"""
function switch_jobs!(inst1, inst2, jobindex1, jobindex2)
    for key in [:g, :n, :o]
        # Swap values between instances
        inst1[key][jobindex1], inst2[key][jobindex2] = inst2[key][jobindex2], inst1[key][jobindex1]
    end
    @debug "Jobs switched between instances"
end


"""
	partition_jobs(n, o, g, p, Pg)

Randomly partitions jobs into `Pg` groups, each satisfying the minimum shelf constraint.

Each group must contain at least `p` subjobs (molds) to meet the shelf availability constraint.
Jobs are assigned randomly to groups until a valid partition is found. This is the initial
partitioning step before the Split-Solve-Merge algorithm begins.

# Arguments
- `n::Vector{Int}`: Job quantities (number of items to produce per job).
- `o::Vector{Int}`: Number of molds per job.
- `g::Vector{Int}`: Job identifiers.
- `p::Int`: Minimum number of subjobs (molds) required per group (number of available shelves).
- `Pg::Int`: Number of groups to create (partition count).

# Returns
- `ns::Tuple{Vector{Int},...}`: Tuple of `Pg` vectors, each containing job quantities for one group.
- `os::Tuple{Vector{Int},...}`: Tuple of `Pg` vectors, each containing mold counts for one group.
- `gs::Tuple{Vector{Int},...}`: Tuple of `Pg` vectors, each containing job IDs for one group.

# Throws
- `error`: If `sum(o) < p × Pg`, partition is infeasible (not enough total molds).

# Algorithm
1. Verify total molds ≥ required molds (`sum(o) ≥ p × Pg`)
2. Randomly assign each job to one of `Pg` groups
3. Check if all groups have ≥ `p` molds
4. Repeat steps 2-3 until valid partition found

# Notes
- Uses random assignment, so results vary between calls
- Used by [`simulated_annealing_init`](@ref) to create initial partition
- See also: [`generate_instance_groups`](@ref), [`simulated_annealing_init`](@ref)

# Examples
```julia
n = [100, 80, 60, 50]
o = [2, 3, 2, 1]  # Total molds: 8
g = [1, 2, 3, 4]
ns, os, gs = partition_jobs(n, o, g, 3, 2)  # 2 groups, each ≥ 3 molds
# Possible result: Group 1 has jobs [1,3], Group 2 has jobs [2,4]
```
"""
function partition_jobs(n, o, g, p, Pg)
    total_subjobs = sum(o)
    required_subjobs = p * Pg
    
    # Check if partition is feasible
    if total_subjobs < required_subjobs
        error("Partition impossible: total subjobs ($total_subjobs) < required ($required_subjobs = $p × $Pg). Need at least $p subjobs per group for $Pg groups: EXIT")
    end
        
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
            @info "Successfully partitioned jobs in $Pg groups after $attempt attempts"
            @debug "Group sizes (subjobs): $sizes"
            return ns, os, gs
        end
    end
end



function generate_instance_groups(instances, Pg)
    if 1 ≤ Pg ≤ 3
        return collect(permutations(instances[1:Pg]))
    else
        error("Partition size not supported: EXIT")
        return []
    end
end


"""
	simulated_annealing_init(; order_dict)

Initializes the simulated annealing process for Split-Solve-Merge optimization.

Creates an initial job partition, generates all permutations of processing stage orders,
solves each permutation using the Split-Solve-Merge algorithm, and returns the best
initial solution as the starting point for simulated annealing.

# Arguments
- `order_dict::Dict`: Configuration parameters with keys:
  - `:g` - Job identifiers
  - `:n` - Job quantities
  - `:o` - Molds per job
  - `:p` - Available shelves
  - `:α`, `:β` - Cost weights
  - `:Nit` - SA iterations
  - `:Pg` - Number of partition groups (stages)

# Returns
- `instances::Vector{Dict}`: Best instance partition (one dict per stage).
- `cost::Float64`: Best total cost among all initial permutations.
- `sols_dict::Vector{Dict}`: Detailed solutions for the best permutation.

# Algorithm
1. Call [`partition_jobs`](@ref) to create initial random partition
2. Generate all permutations of stages via [`generate_instance_groups`](@ref)
3. Solve each permutation with [`process_instances_group`](@ref)
4. Return permutation with minimum total cost

# Notes
- For `Pg = 1`: Only one instance, no permutations (exact MILP solve)
- For `Pg = 2`: 2 permutations (A→B, B→A)
- For `Pg = 3`: 6 permutations (all orderings of 3 stages)
- Logs all permutation costs via `@info` for debugging
- See also: [`partition_jobs`](@ref), [`generate_instance_groups`](@ref), [`simulated_annealing`](@ref)

# Examples
```julia
order_dict = Dict(:g => [1,2,3], :n => [100,80,60], :o => [2,3,2], :p => 4, :Pg => 2, ...)
instances, cost, sols = simulated_annealing_init(; order_dict=order_dict)
println("Initial cost: ", cost)
```
"""
function simulated_annealing_init(; order_dict)
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

    costs, sols_dict_group = process_instances_group(instances_group; order_dict=order_dict)

    #@info "Sols_dict_group $sols_dict_group" 
    @info "Costs  $costs"
    cost, min_indx = findmin(costs)
    @info "Best cost $cost at instances $min_indx"
    instances = instances_group[min_indx]
    @debug "Instances  $instances"
    sols_dict = sols_dict_group[min_indx]

    for (i, sols) ∈ enumerate(sols_dict_group)
        @debug "Instance permutation $i"
        @debug sols
    end

    return instances, cost, sols_dict

end


"""
	neighbour_partition!(instances, g, p, Pg)

Generates a neighboring job partition by moving or swapping jobs between stages.

This function implements the neighborhood structure for simulated annealing optimization.
It modifies the current partition in-place by attempting to move or swap jobs while
preserving the constraint that each stage must have at least `p` molds.

# Arguments
- `instances::Vector{Dict}`: Current partition (one dict per stage), modified in-place.
- `g::Vector{Int}`: Set of all possible job IDs in the problem.
- `p::Int`: Minimum molds required per stage (shelf constraint).
- `Pg::Int`: Number of partition groups (must be 2 or 3).

# Returns
- `jobid::Int`: The job ID that was moved or swapped.

# Throws
- `error`: If no valid neighbor found after `length(g)` attempts (partition is stuck).
- `error`: If `Pg < 2` or `Pg > length(instances)` (unsupported partition size).

# Algorithm
For each job ID sampled randomly from `g`:
1. Try all permutations of stages as (source, destination) pairs
2. **Move Strategy** (if `sum(src[:o]) > p`):
   - If removing job leaves ≥ `p` molds in source, move job to destination
   - Otherwise, rebalance by moving jobs from destination to source until feasible
3. **Swap Strategy** (if `sum(src[:o]) == p`):
   - Find matching job in destination with same mold count
   - Swap jobs to maintain feasibility in both stages

# Notes
- The `!` suffix indicates in-place modification of `instances`
- Maximum attempts: `length(g)` (one try per job)
- Uses [`move_job!`](@ref) and [`switch_jobs!`](@ref) for partition modifications
- See also: [`simulated_annealing`](@ref), [`move_job!`](@ref), [`switch_jobs!`](@ref)

# Examples
```julia
instances = [
    Dict(:g => [1,2], :o => [2,3], :n => [100,80]),
    Dict(:g => [3,4], :o => [2,2], :n => [60,50])
]
jobid = neighbour_partition!(instances, [1,2,3,4], 4, 2)
# instances modified: job moved/swapped between stages
```
"""
function neighbour_partition!(instances, g, p, Pg)

    if Pg < 2 || Pg > length(instances)
        error("Partition size not supported: EXIT")
    end

    counter = 0
    max_counter = length(g)
    #=
    		while counter < max_counter
    			counter += 1
    			@debug "Counter $counter for neighbour_partition"
    			jobid = rand(g)
    			for partition ∈ collect(permutations(instances[1:Pg]))
    				src = partition[1]
    				if length(partition) < 2
    					@debug "Only one partition available, skipping"
    					continue
    				end
    				dest = rand(partition[2:end])  # Randomly pick one of the remaining destinations
    				# Warning: Some restrictions to avoid moving/swaping jobs when multiple molds are available
    				if jobid ∈ src[:g]
    					if sum(src[:o]) > p
    						jobsrcindex = findfirst(==(jobid), src[:g])
    						if jobsrcindex === nothing
    							continue  # Try another partition
    						end
    						neighdiff = sum(src[:o]) - src[:o][jobsrcindex]
    						if neighdiff ≥ p
    							move_job!(src, dest, jobsrcindex)
    							return jobid
    						elseif neighdiff < p
    							countdiff = p - neighdiff
    							while countdiff > 0
    								eligible_jobs = findall(x -> dest[:o][x] ≤ countdiff, 1:length(dest[:o]))
    								if !isempty(eligible_jobs)
    									chosen_idx = rand(eligible_jobs)
    									move_job!(dest, src, chosen_idx)
    									countdiff -= dest[:o][chosen_idx]
    									if countdiff == 0
    										return jobid
    									end
    								else
    									break
    								end
    							end
    						end
    					elseif sum(src[:o]) == p
    						jobsrcindex = findfirst(==(jobid), src[:g])
    						if jobsrcindex === nothing
    							continue  # Try another partition
    						end
    						# Find all jobs in dest such that dest[:o][jobdestindex] == src[:o][jobsrcindex]
    						matching_dest_indices = findall(x -> dest[:o][x] == src[:o][jobsrcindex], 1:length(dest[:o]))
    						if !isempty(matching_dest_indices)
    							# Choose randomly one of them
    							jobdestindex = rand(matching_dest_indices)
    							# Switch jobs
    							switch_jobs!(src, dest, jobsrcindex, jobdestindex)
    							return jobid
    						end
    					end
    				end
    			end
    			=#

    while counter < max_counter
        counter += 1
        @debug "Counter $counter for neighbour_partition"

        jobid = rand(g)

        for partition ∈ collect(permutations(instances[1:Pg]))
            if length(partition) < 2
                continue  # not enough instances to consider a source and a destination
            end

            src = partition[1]
            dest_candidates = partition[2:end]
            if isempty(dest_candidates)
                continue
            end

            dest = rand(dest_candidates)

            # Check that src has jobid
            jobsrcindex = findfirst(==(jobid), src[:g])
            if jobsrcindex === nothing
                continue
            end

            # CASE 1: src has more than p subjobs → direct move
            if sum(src[:o]) > p
                neighdiff = sum(src[:o]) - src[:o][jobsrcindex]

                if neighdiff ≥ p
                    move_job!(src, dest, jobsrcindex)
                    @debug "Moved job $jobid from src to dest"
                    return jobid
                else
                    # Try to rebalance from dest to src
                    countdiff = p - neighdiff
                    while countdiff > 0
                        eligible_jobs = findall(x -> dest[:o][x] ≤ countdiff, 1:length(dest[:o]))
                        if isempty(eligible_jobs)
                            break
                        end

                        chosen_idx = rand(eligible_jobs)
                        if chosen_idx > length(dest[:o])
                            @debug "chosen_idx out of bounds: $chosen_idx for dest[:o] of length $(length(dest[:o]))"
                            break
                        end
                        countdiff -= dest[:o][chosen_idx]
                        move_job!(dest, src, chosen_idx)
                        if countdiff ≤ 0
                            @debug "Rebalanced and moved job $jobid"
                            return jobid
                        end
                    end
                end

                # CASE 2: src has exactly p subjobs → try swap
            elseif sum(src[:o]) == p
                # Ensure jobsrcindex is valid
                if jobsrcindex > length(src[:o])
                    @debug "jobsrcindex $jobsrcindex out of bounds for src[:o] of length $(length(src[:o]))"
                    continue
                end

                src_job_o = src[:o][jobsrcindex]
                matching_dest_indices = findall(x -> dest[:o][x] == src_job_o, 1:length(dest[:o]))

                if isempty(matching_dest_indices)
                    continue
                end

                jobdestindex = rand(matching_dest_indices)

                if jobdestindex > length(dest[:g])
                    @debug "jobdestindex $jobdestindex out of bounds for dest[:g] of length $(length(dest[:g]))"
                    continue
                end

                switch_jobs!(src, dest, jobsrcindex, jobdestindex)
                @debug "Switched job $jobid with dest job at index $jobdestindex"
                return jobid
            end
        end

        # If after all permutations no move is possible

        if counter == max_counter
            error("No available neighbours $max_counter attempts: EXIT")
        end
    end
end


"""
	simulated_annealing(input_instances, input_cost, input_sols_dict; order_dict)

Performs simulated annealing to optimize job partitioning in the Split-Solve-Merge algorithm.

This function iteratively explores the solution space by generating neighboring partitions,
evaluating their costs via the Split-Solve-Merge algorithm, and accepting or rejecting them
based on the Metropolis criterion with exponential cooling.

# Arguments
- `input_instances::Vector{Dict}`: Initial partition from [`simulated_annealing_init`](@ref).
- `input_cost::Float64`: Initial total cost.
- `input_sols_dict::Vector{Dict}`: Initial detailed solutions.
- `order_dict::Dict`: Configuration with keys:
  - `:g`, `:p`, `:Pg` - Job parameters
  - `:Nit` - Number of SA iterations
  - `:T0`, `:Tf`, `:Tj` - Temperature parameters (initial, final, jump frequency)
  - `:Gl` - Global time limit in seconds

# Returns
- `best_instances::Vector{Dict}`: Best partition found during search.
- `best_cost::Float64`: Best total cost achieved.
- `best_sols_dict::Vector{Dict}`: Detailed solutions for best partition.
- `cur_costs_list::Vector{Float64}`: Cost trajectory over iterations (for convergence analysis).

# Algorithm
1. **Cooling Schedule**: Exponential cooling `T = T₀ × γⁱ` where `γ = (Tf/T₀)^(1/Nit)`, updated every `Tj` iterations
2. **Iteration Loop**:
   - Generate neighbor via [`neighbour_partition!`](@ref)
   - Evaluate all permutations via [`process_instances_group`](@ref)
   - Compute Δ = (cur_cost - alt_cost) / T
   - **Accept** if alt_cost ≤ cur_cost
   - **Metropolis Accept** if exp(Δ) > rand() (uphill moves)
   - Update best if improvement found
3. **Termination**: After `Nit` iterations or `Gl` seconds elapsed

# Notes
- Uses deep copies to preserve solution history
- Logs detailed iteration info: costs, temperature, acceptance probability
- Tracks elapsed time and stops early if global time limit exceeded
- See also: [`simulated_annealing_init`](@ref), [`neighbour_partition!`](@ref), [`run`](@ref)

# Examples
```julia
init_inst, init_cost, init_sols = simulated_annealing_init(; order_dict=od)
best_inst, best_cost, best_sols, costs = simulated_annealing(
    init_inst, init_cost, init_sols; order_dict=od
)
println("Improvement: ", init_cost - best_cost)
```
"""
function simulated_annealing(input_instances, input_cost, input_sols_dict; order_dict)
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
        alt_costs, alt_sols_dict_group = process_instances_group(instances_group; order_dict=order_dict)

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
	run(order_dict)

Main entry point for the Split-Solve-Merge algorithm with simulated annealing optimization.

Executes the complete Split-Solve-Merge workflow: initializes job partitioning, runs simulated
annealing to optimize partition ordering (if `Pg > 1`), and logs all results to disk with
appropriate directory structure based on solution method (exact vs. heuristic).

# Arguments
- `order_dict::Dict`: Complete configuration dictionary with keys:
  - `:g`, `:n`, `:o` - Job data (IDs, quantities, molds)
  - `:p`, `:α`, `:β` - Problem parameters
  - `:Pg` - Partition groups (1 = exact MILP, >1 = heuristic Split-Solve-Merge)
  - `:Nit`, `:T0`, `:Tf`, `:Tj`, `:Gl` - SA parameters
  - `:Oid` - Order identifier for file naming

# Returns
- `cost::Float64`: Final optimized total cost.
- `m::Int`: Total number of machines used across all stages.
- `elapsed_time::Float64`: Total execution time in seconds.

# Algorithm
1. **Initialization**: Call [`simulated_annealing_init`](@ref) to create initial partition and evaluate permutations
2. **Single Stage (Pg = 1)**: Return immediately (exact MILP, no SA needed)
3. **Multi-Stage (Pg > 1)**: Run [`simulated_annealing`](@ref) to optimize partition
4. **Cost Aggregation**: Sum costs across stages, accounting for job transfers (use `:cost_p` when jobs transferred)
5. **Logging**: Save all results to `data/sims/exact/` or `data/sims/heuristics/` via demux logger

# Side Effects
- Creates log files in `datadir("sims/heuristics", Oid)` or `datadir("sims/exact", Oid)`
- Logs include: order parameters, partition details, solution matrices, costs, timings

# Notes
- Uses DrWatson's `datadir()` for reproducible directory structure
- Logs via `demux_logger()` for multiple output streams
- Machine count `m` calculated by checking shelf utilization in last machines
- See also: [`simulated_annealing_init`](@ref), [`simulated_annealing`](@ref)

# Examples
```julia
order_dict = Dict(
    :g => [1,2,3], :n => [100,80,60], :o => [2,3,2],
    :p => 4, :α => 1.0, :β => 10.0, :Pg => 2, :Nit => 100,
    :T0 => 100.0, :Tf => 0.01, :Tj => 5, :Gl => 3600.0,
    :Oid => "order_001"
)
cost, machines, time = run(order_dict)
println("Total cost: \$cost using \$machines machines in \$time seconds")
```
"""
function run(order_dict)
    prefix = nothing
    if order_dict[:Pg] > 1
        prefix = datadir("sims/heuristics", order_dict[:Oid])
    else
        prefix = datadir("sims/exact",  order_dict[:Oid])
    end

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
            return init_cost, init_sols_dict[1][:m], elapsed_time
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
                m = round(Int, sol[:m])
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
        return best_cost, h_m, elapsed_time
    end
end

# End of module
end
