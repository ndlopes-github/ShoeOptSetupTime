#= Copyright (C) 2024
Nuno David Lopes.
Created:  2024/04/09
Last changed - N. Lopes:2025/07/26 16:26:19
=#
using DrWatson
@quickactivate "ShoeOptSetupTime"

# Load settings FIRST, before including the module
include(datadir("settings", "E_O2_#2_3p.jl"))

# Now include and use the module
include(scriptsdir("split_solve_merge_milp.jl"))
using .SplitSolveMergeMILP

# Pass the already-constructed dictionary
SplitSolveMergeMILP.run(order_dict)

