#= Copyright (C) 2024
Nuno David Lopes.
Created:  2026/04/20
Last changed - N. Lopes:2026/04/20

CLI runner for SSM-SA (Split-Solve-Merge Simulated Annealing).
Loads a settings instance, overrides SA parameters with irace Config 4 tuned values
(T0=4.5479, Tf=0.0373, Nit=150), and runs SplitSolveMergeMILP.run().

Usage:
    julia --project scripts/run_ssm_sa_cli.jl --instance=H_O1_#1_2p.jl
    julia --project scripts/run_ssm_sa_cli.jl --instance H_O1_#1_2p.jl
=#
using DrWatson
@quickactivate "SoftIdea"

# ─── Irace Config 4 tuned parameters ─────────────────────────────────────────
const IRACE_T0 = 4.5479
const IRACE_Tf = 0.0373
const IRACE_Nit = 150

"""
Parse CLI arguments and return the instance settings filename.

# Arguments
- `args`: Command-line arguments vector

# Returns
- `String`: Instance settings filename
"""
function parse_instance_arg(args::Vector{String})::String
    positional_args = String[]

    for (index, arg) in pairs(args)
        if arg == "--instance"
            index == lastindex(args) && error("Missing value for --instance: EXIT")
            return strip(args[index+1])
        end

        if startswith(arg, "--instance=")
            value = strip(split(arg, "=", limit=2)[2])
            isempty(value) && error("Missing value for --instance: EXIT")
            return value
        end

        startswith(arg, "--") || push!(positional_args, strip(arg))
    end

    isempty(positional_args) && error(
        "No instance argument provided. Use --instance <file>, --instance=<file>, or a positional file name: EXIT",
    )

    return first(positional_args)
end

"""Parse optional --beta argument. Returns the value as a Float64, or nothing."""
function parse_beta_arg(args::Vector{String})
    for (index, arg) in pairs(args)
        if arg == "--beta"
            index == lastindex(args) && error("Missing value for --beta: EXIT")
            return parse(Float64, args[index+1])
        end
        if startswith(arg, "--beta=")
            return parse(Float64, split(arg, "=", limit=2)[2])
        end
    end
    return nothing
end

instance = parse_instance_arg(ARGS)
beta_override = parse_beta_arg(ARGS)
settings_path = datadir("settings", instance)
isfile(settings_path) || error("Settings file not found: $(settings_path): EXIT")

@info "Running SSM-SA for instance: $instance"
@info "Parameters (irace Config 4): T0=$(IRACE_T0)  Tf=$(IRACE_Tf)  Nit=$(IRACE_Nit)"
!isnothing(beta_override) && @info "Beta override: β = $beta_override (instance value bypassed)"

include(settings_path)
# order_dict is now defined in the current scope via include

# ─── Override with irace-tuned parameters ────────────────────────────────────
order_dict[:T0] = IRACE_T0
order_dict[:Tf] = IRACE_Tf
order_dict[:Nit] = IRACE_Nit
!isnothing(beta_override) && (order_dict[:β] = beta_override)

# Update Oid to reflect the tuned parameters
const FILEBASENAME = splitext(basename(settings_path))[1]
order_dict[:Oid] = "$(FILEBASENAME)_irace_T0_$(IRACE_T0)_Tf_$(IRACE_Tf)_Nit_$(IRACE_Nit)"

include(scriptsdir("split_solve_merge_milp.jl"))
using .SplitSolveMergeMILP

SplitSolveMergeMILP.run(order_dict)
