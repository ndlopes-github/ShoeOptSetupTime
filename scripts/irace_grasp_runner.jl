#= Copyright (C) 2024
Nuno David Lopes.
Irace target runner for GRASP.

Called as:
    julia --project=<project> scripts/irace_grasp_runner.jl \
        <instance_file> <seed> --Nit <v>

Runs GRASP for Nit independent constructions and prints the best cost.
All logging is suppressed so that only the cost value appears on stdout.
=#

using DrWatson
@quickactivate "SoftIdea"

using Logging, Random

include(scriptsdir("simulated_annealing.jl"))
using .SimulatedAnnealing: get_cost_from_shelves

include(scriptsdir("grasp.jl"))
using .Grasp

# ─── Parse arguments ─────────────────────────────────────────────────────────
if length(ARGS) < 2
    error("Usage: irace_grasp_runner.jl <instance_file> <seed> [--param value ...]: EXIT")
end

instance_file = abspath(ARGS[1])
seed = parse(Int, ARGS[2])

param_args = ARGS[3:end]
function get_param(args, flag)
    idx = findfirst(==(flag), args)
    isnothing(idx) ? nothing : args[idx+1]
end

Nit_arg = get_param(param_args, "--Nit")

# ─── Load instance ────────────────────────────────────────────────────────────
if !isfile(instance_file)
    error("Instance file not found: $(instance_file): EXIT")
end

with_logger(NullLogger()) do
    include(instance_file)
end

# order_dict is now defined in the current scope via include
if !isnothing(Nit_arg)
    order_dict[:Nit] = parse(Int, Nit_arg)
end

# ─── Seed and run ─────────────────────────────────────────────────────────────
Random.seed!(seed)

@unpack α, β = order_dict

best_shelves = redirect_stdout(devnull) do
    with_logger(NullLogger()) do
        Grasp.run(order_dict; log_every=0)
    end
end
_, _, best_cost = get_cost_from_shelves(best_shelves, α, β)

# ─── Output ───────────────────────────────────────────────────────────────────
# Irace expects exactly one number on stdout
println(best_cost)
