#= Copyright (C) 2024
Nuno David Lopes.
Irace target runner for pure Simulated Annealing.

Called as:
    julia --project=<project> scripts/irace_sa_runner.jl \
        <instance_file> <seed> --T0 <v> --Tf <v> --Nit <v> --Tj <v> --Nruns <v>

Prints a single line: the minimum cost achieved across Nruns independent SA runs.
All SA logging is suppressed so that only the cost value appears on stdout.
=#

using DrWatson
@quickactivate "SoftIdea"

using Logging, Random

include(scriptsdir("simulated_annealing.jl"))
using .SimulatedAnnealing

# ─── Parse arguments ─────────────────────────────────────────────────────────
if length(ARGS) < 2
    error("Usage: irace_sa_runner.jl <instance_file> <seed> [--param value ...]: EXIT")
end

instance_file = abspath(ARGS[1])
seed = parse(Int, ARGS[2])

# Parse optional named parameters
param_args = ARGS[3:end]
function get_param(args, flag)
    idx = findfirst(==(flag), args)
    isnothing(idx) ? nothing : args[idx+1]
end

T0_arg = get_param(param_args, "--T0")
Tf_arg = get_param(param_args, "--Tf")
Nit_arg = get_param(param_args, "--Nit")
Tj_arg = get_param(param_args, "--Tj")
Nruns_arg = get_param(param_args, "--Nruns")

# ─── Load instance ────────────────────────────────────────────────────────────
if !isfile(instance_file)
    error("Instance file not found: $(instance_file): EXIT")
end

# Suppress @info/@debug output from instance file and SA
with_logger(NullLogger()) do
    include(instance_file)
end

# order_dict is now defined in the current scope via include
# Override SA parameters if provided on the command line
if !isnothing(T0_arg)
    order_dict[:T0] = parse(Float64, T0_arg)
end
if !isnothing(Tf_arg)
    order_dict[:Tf] = parse(Float64, Tf_arg)
end
if !isnothing(Nit_arg)
    order_dict[:Nit] = parse(Int, Nit_arg)
end
if !isnothing(Tj_arg)
    order_dict[:Tj] = parse(Int, Tj_arg)
end

Nruns = isnothing(Nruns_arg) ? 100 : parse(Int, Nruns_arg)

# ─── Seed and run ─────────────────────────────────────────────────────────────
Random.seed!(seed)

@unpack α, β = order_dict

function _run_trials(order_dict, α, β, nruns)
    best = Inf
    for _ in 1:nruns
        shelves = SimulatedAnnealing.run(order_dict; log_every=0)
        _, _, cost = SimulatedAnnealing.get_cost_from_shelves(shelves, α, β)
        cost < best && (best = cost)
    end
    return best
end

best_cost = with_logger(NullLogger()) do
    _run_trials(order_dict, α, β, Nruns)
end

# ─── Output ───────────────────────────────────────────────────────────────────
# Irace expects exactly one number on stdout
println(best_cost)
