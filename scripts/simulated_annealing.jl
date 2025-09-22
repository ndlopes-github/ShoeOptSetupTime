#= Copyright (C) 2024
Nuno David Lopes.
Created:  2024/04/09
Last changed - N. Lopes:2025/07/26 16:32:21
=#
module SimulatedAnnealing
export run_sim

using DrWatson
@quickactivate "SoftIdea"
using Logging, Dates, DataFrames, Gurobi, Random, JuMP, Combinatorics, Printf

include(srcdir("loggers.jl"))



"""
	model_builder(instance; x_cnst_indxs = []; order_dict = order_dict)

Builds the optimization model using the given instance and optional constraints.

# Arguments
- `instance`: The instance data for the model.
- `x_cnst_indxs`: Optional constraints for the model.

# Returns
- `model`: The built optimization model.
"""
function model_builder(instance; x_cnst_indxs = [], order_dict = order_dict)
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
	@variable(model, t[i ∈ Is, j ∈ Js, k ∈ Ks[j]] ≤ sum(n), Int)
	@variable(model, tbar[i ∈ Is])
	@variable(model, y[i ∈ Is])

	@objective(model, Min, sum(α * tbar[i] + β * y[i] for i in Is))
	@constraint(model, [j ∈ Js], sum(t[i, j, k] for i ∈ Is, k ∈ Ks[j]) == n[j])
	@constraint(model, [i ∈ Is, j ∈ Js, k ∈ Ks[j]], tbar[i] ≥ t[i, j, k])
	@constraint(model, [i in Is], sum(x[i, j, k] for j in Js, k ∈ Ks[j]) ≤ p)
	@constraint(model, [i ∈ Is, j ∈ Js, k ∈ Ks[j]], t[i, j, k] ≤ n[j] * x[i, j, k])
	@constraint(model, [i ∈ Is, j ∈ Js, k ∈ Ks[j]], t[i, j, k] ≥ x[i, j, k])

	#@constraint(model,
	#    [i ∈ Is, j ∈ Js, l ∈ Js, u ∈ Ks[j], v ∈ Ks[l]],
	#    t[i, j, u] ≤ t[i, l, v] + n[j] * (1 - x[i, l, v]))

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
	solver(model, instance; order_dict = order_dict)

Solves the given optimization model and returns the solution.

# Arguments
- `model`: The optimization model to solve.
- `instance`: The instance data for the model.

# Returns
- `table_dict`: A dictionary with the solution summary.
- `sol_dict`: A dictionary with the detailed solution.
"""
function solver(model, instance; order_dict = order_dict)
	@debug "Solver for Instance:  $instance "
	@unpack n, o, g = instance
	@unpack α, β, p = order_dict

	optimize!(model)
	MiB = round(Sys.maxrss() / 1048576; digits = 1)


	if !is_solved_and_feasible(model; dual = false)
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

	stime = round(solve_time(model); digits = 10)
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
	determine_transfer_indexes(df, sols_dict, i2)

Determines the indexes of jobs to transfer based on the solution.

# Arguments
- `df`: The DataFrame with the solution summary.
- `sols_dict`: The dictionary with the detailed solution.
- `i2`: The instance to transfer jobs to.

# Returns
- `indexes`: The indexes of jobs to transfer.
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
	solve_and_store(instance, df, sols_dict, x_cnst_indxs = []; order_dict = order_dict)

Solves the given instance and stores the solution in the provided DataFrame and dictionary.

# Arguments
- `instance`: The instance data for the model.
- `df`: The DataFrame to store the solution summary.
- `sols_dict`: The dictionary to store the detailed solution.
- `x_cnst_indxs`: Optional constraints for the model.

# Returns
- `df`: The updated DataFrame with the solution summary.
- `sols_dict`: The updated dictionary with the detailed solution.
"""
function solve_and_store(instance, df, sols_dict, x_cnst_indxs = []; order_dict = order_dict)
	model = model_builder(instance; x_cnst_indxs = x_cnst_indxs, order_dict = order_dict)
	table_dict, sol_dict = solver(model, instance; order_dict = order_dict)

	if table_dict !== nothing
		push!(df, only(DataFrame(table_dict)))
		push!(sols_dict, sol_dict)
	end

	@debug " Instance solved and stored"
	return df, sols_dict
end



"""
	process_instances_group(instances_group)

Processes a group of instances and returns the costs and solutions.

# Arguments
- `instances_group`: The group of instances to process.

# Returns
- `costs`: The list of costs for each instance.
- `sols_dict_group`: The list of solution dictionaries for each instance.
"""
function process_instances_group(instances_group; order_dict = order_dict)
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
		df, sols_dict = solve_and_store(inst[1], df, sols_dict; order_dict = order_dict)

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
			df, sols_dict = solve_and_store(inst[j], df, sols_dict, x_indxs; order_dict = order_dict)
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

Transfers reminder jobs, if any, to the destination instance.

# Arguments
- `dest`: The destination instance.
- `indexes`: The indexes of jobs to transfer.
- `val`: The value to assign to the transferred jobs.

# Returns
- `x_indxs`: The updated list of indexes.
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

Moves a job from the source instance to the destination instance.

# Arguments
- `src`: The source instance.
- `dest`: The destination instance.
- `jobindex`: The index of the job to move.
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

Swaps the jobs at `jobindex1` in `inst1` with `jobindex2` in `inst2`.

# Arguments
- `inst1`: The first instance.
- `inst2`: The second instance.
- `jobindex1`: The index of the job in `inst1` to be swapped.
- `jobindex2`: The index of the job in `inst2` to be swapped.
"""
function switch_jobs!(inst1, inst2, jobindex1, jobindex2)
	for key in [:g, :n, :o]
		# Swap values between instances
		inst1[key][jobindex1], inst2[key][jobindex2] = inst2[key][jobindex2], inst1[key][jobindex1]
	end
	@debug "Jobs switched between instances"
end


"""
	partition_jobs(n, o, g, p)

Partitions the jobs into two groups such that each group has at least `p` subjobs.

# Arguments
- `n`: A list of job sizes.
- `o`: A list of subjob sizes corresponding to each job.
- `g`: A list of job identifiers.
- `p`: The minimum number of subjobs required in each group. (number of available shelves)

# Returns
- `ns`: A tuple containing two lists of job sizes, one for each group.
- `os`: A tuple containing two lists of subjob sizes, one for each group.
- `gs`: A tuple containing two lists of job identifiers, one for each group.
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



function generate_instance_groups(instances, Pg)
	if 1 ≤ Pg ≤ 3
		return collect(permutations(instances[1:Pg]))
	else
		error("Partition size not supported: EXIT")
		return []
	end
end


"""
	simulated_annealing_init(order_dict)

Initializes the simulated annealing process.

# Arguments
- `order_dict`: The dictionary with the order data.

# Returns
- `instances`: The initialized instances.
- `cost`: The initial cost.
- `sols_dict`: The initial solutions dictionary.
"""
function simulated_annealing_init(; order_dict = order_dict)
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

	for (i, sols) ∈ enumerate(sols_dict_group)
		@debug "Instance permutation $i"
		@debug sols
	end

	return instances, cost, sols_dict

end


"""
Finds a valid job partition and moves it.

# Arguments
- `instances`: The current job instances.
- `g`: Set of possible job IDs.
- `p`: Partition size limit.

# Returns
- `jobid`: The selected job ID.
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
	simulated_annealing!(best_instances, best_cost, best_sols_dict, order_dict)

Performs an iteration of the simulated annealing process.

# Arguments
- `best_instances`: The best instances so far.
- `best_cost`: The best cost so far.
- `best_sols_dict`: The best solutions dictionary so far.
- `order_dict`: The dictionary with the order data.

# Returns
- `cur_costs_list`: The list of current costs after the iteration.
"""
function simulated_annealing(input_instances, input_cost, input_sols_dict; order_dict = order_dict)
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
	run_sim(;order_file = nothing)

Runs the simulated annealing process for the given order dictionary and logs the results.

# Arguments
- `order_file`: The order file to use.
"""
function run_sim(; order_file = nothing)
	if order_file == nothing
		error("No order file provided: EXIT")
		return nothing
	end

	include(datadir("settings", order_file))
	prefix = nothing
	if order_dict[:Pg] > 1
		prefix = datadir("sims/heuristics", order_file[1:end-3] * order_dict[:Oid])
	else
		prefix = datadir("sims/exact", order_file[1:end-3] * order_dict[:Oid])
	end

	with_logger(demux_logger(prefix)) do
		elapsed_time = @elapsed begin
			@info "Simulation date  $(now())"
			@info "Order Dictionary $(order_dict)"

			init_instances, init_cost, init_sols_dict = simulated_annealing_init(; order_dict = order_dict)

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
			best_instances, best_cost, best_sols_dict, cur_costs_list = simulated_annealing(init_instances, init_cost, init_sols_dict; order_dict = order_dict)
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

