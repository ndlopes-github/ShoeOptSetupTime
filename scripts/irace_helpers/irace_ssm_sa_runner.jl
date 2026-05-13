#= Copyright (C) 2024
Nuno David Lopes.
irace target runner for Split-Solve-Merge Simulated Annealing (SSM-SA).

Called by irace as:
    julia ... scripts/irace_ssm_sa_runner.jl \
        <instance_file> <seed> <candidate_id> --T0 <v> --Tf <v> --Nit <v>

Prints a single line: the objective cost returned by SplitSolveMergeMILP.run().
All solver/SA logging is suppressed so only the cost appears on stdout.

Notes:
  - No Nruns: SSM-SA runs a single trajectory per call (MILP overhead precludes multi-run)
  - No Tj: geometric cooling schedule is determined by T0, Tf, Nit internally
  - Pg, Tl, Gl are fixed at their paper values (Pg=2, Tl=30, Gl=1800)
  - Oid is made unique per call using candidate_id + seed to avoid log file conflicts
    during parallel irace execution
  - MILP-limited instances (e.g. H_O2_#3_2p) will produce a valid cost even when the
    MILP sub-solver hits its time limit — irace receives the cost and can evaluate it
=#

using DrWatson
@quickactivate "SoftIdea"

using Logging, Random

include(scriptsdir("split_solve_merge_milp.jl"))
using .SplitSolveMergeMILP

# ─── Parse arguments ─────────────────────────────────────────────────────────
if length(ARGS) < 3
    error("Usage: irace_ssm_sa_runner.jl <instance_file> <seed> <candidate_id> [--param value ...]: EXIT")
end

instance_file  = abspath(ARGS[1])
seed           = parse(Int, ARGS[2])
candidate_id   = ARGS[3]

param_args = ARGS[4:end]

function get_param(args, flag)
    idx = findfirst(==(flag), args)
    isnothing(idx) ? nothing : args[idx+1]
end

T0_arg  = get_param(param_args, "--T0")
Tf_arg  = get_param(param_args, "--Tf")
Nit_arg = get_param(param_args, "--Nit")

# ─── Load instance ────────────────────────────────────────────────────────────
if !isfile(instance_file)
    error("Instance file not found: $(instance_file): EXIT")
end

with_logger(NullLogger()) do
    include(instance_file)
end
# order_dict is now defined in the current scope via include

# ─── Override tuned parameters ───────────────────────────────────────────────
if !isnothing(T0_arg)
    order_dict[:T0]  = parse(Float64, T0_arg)
end
if !isnothing(Tf_arg)
    order_dict[:Tf]  = parse(Float64, Tf_arg)
end
if !isnothing(Nit_arg)
    order_dict[:Nit] = parse(Int, Nit_arg)
end

# Make Oid unique per irace call to prevent log file collisions during parallel runs
order_dict[:Oid] = "irace_ssm_sa_cid$(candidate_id)_seed$(seed)_T0$(order_dict[:T0])_Tf$(order_dict[:Tf])_Nit$(order_dict[:Nit])"

# ─── Seed and run ─────────────────────────────────────────────────────────────
Random.seed!(seed)

# Redirect C-level stdout (Gurobi, MILP solver) to devnull so only the cost
# reaches irace. Julia's own logging is already suppressed via NullLogger.
cost = redirect_stdout(devnull) do
    with_logger(NullLogger()) do
        best_cost, _m, _elapsed = SplitSolveMergeMILP.run(order_dict)
        best_cost
    end
end

# ─── Output ───────────────────────────────────────────────────────────────────
# irace expects exactly one number on stdout
println(cost)
