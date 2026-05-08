#= Copyright (C) 2024
Nuno David Lopes.
Created:  2026/03/19
Last changed - N. Lopes:2026/03/19 17:46:24
=#
using DrWatson
@quickactivate "SoftIdea"

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
            return strip(args[index + 1])
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

instance = parse_instance_arg(ARGS)
settings_path = datadir("settings", instance)
isfile(settings_path) || error("Settings file not found: $(settings_path): EXIT")

@info "Running MILP for instance: $instance"

include(settings_path)
include(scriptsdir("split_solve_merge_milp.jl"))
using .SplitSolveMergeMILP

SplitSolveMergeMILP.run(order_dict)
