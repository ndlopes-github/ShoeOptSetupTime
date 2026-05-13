#= Copyright (C) 2025
Nuno David Lopes.
Irace target runner for Genetic Algorithm.

Called as:
    julia --project=<project> scripts/irace_ga_runner.jl \
        <instance_file> <seed> --Nit <v> --pop_size <v> \
        --clone_threshold <v> --Nruns <v>

Runs GA for Nruns independent runs (each with Nit generations) and prints
the minimum cost achieved. All logging is suppressed so that only the cost
value appears on stdout.
=#

using DrWatson
@quickactivate "SoftIdea"

using Logging, Random

include(scriptsdir("simulated_annealing.jl"))
using .SimulatedAnnealing: get_cost_from_shelves

include(scriptsdir("genetic_algorithm.jl"))
using .GeneticAlgorithm

# ─── Parse arguments ──────────────────────────────────────────────────────────
if length(ARGS) < 2
    error("Usage: irace_ga_runner.jl <instance_file> <seed> [--param value ...]: EXIT")
end

instance_file = abspath(ARGS[1])
seed = parse(Int, ARGS[2])

param_args = ARGS[3:end]
function get_param(args, flag)
    idx = findfirst(==(flag), args)
    isnothing(idx) ? nothing : args[idx+1]
end

Nit_arg = get_param(param_args, "--Nit")
pop_size_arg = get_param(param_args, "--pop_size")
clone_threshold_arg = get_param(param_args, "--clone_threshold")
Nruns_arg = get_param(param_args, "--Nruns")

# ─── Load instance ────────────────────────────────────────────────────────────
isfile(instance_file) || error("Instance file not found: $(instance_file): EXIT")

with_logger(NullLogger()) do
    include(instance_file)
end

# order_dict is now defined in the current scope via include

# ─── Override GA parameters ───────────────────────────────────────────────────
isnothing(Nit_arg) && error("Missing required argument --Nit: EXIT")
isnothing(pop_size_arg) && error("Missing required argument --pop_size: EXIT")
isnothing(clone_threshold_arg) && error("Missing required argument --clone_threshold: EXIT")
isnothing(Nruns_arg) && error("Missing required argument --Nruns: EXIT")

order_dict[:Nit] = parse(Int, Nit_arg)
pop_size = parse(Int, pop_size_arg)
clone_threshold = parse(Float64, clone_threshold_arg)
Nruns = parse(Int, Nruns_arg)

# ─── Seed and run ─────────────────────────────────────────────────────────────
Random.seed!(seed)

@unpack α, β = order_dict

function _run_ga_trials(order_dict, α, β, nruns, pop_size, clone_threshold)
    best = Inf
    for _ in 1:nruns
        shelves = with_logger(NullLogger()) do
            GeneticAlgorithm.run_ga(order_dict;
                pop_size=pop_size,
                clone_threshold=clone_threshold,
                log_every=0)
        end
        _, _, cost = get_cost_from_shelves(shelves, α, β)
        cost < best && (best = cost)
    end
    return best
end

best_cost = redirect_stdout(devnull) do
    _run_ga_trials(order_dict, α, β, Nruns, pop_size, clone_threshold)
end

# ─── Output ───────────────────────────────────────────────────────────────────
# Irace expects exactly one number on stdout
println(best_cost)
