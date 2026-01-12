#= Copyright (C) 2024
Nuno David Lopes.
Created:  2025/11/21
Last changed - N. Lopes: 2026/01/08 12:08:31
=#

using DrWatson
@quickactivate "ShoeOptSetupTime"

# Load settings FIRST, before including the module
include(datadir("settings", "H_O2_#2_3p.jl"))

# Now include and use the GRASP module's greedy wrapper
include(scriptsdir("grasp.jl"))
using .Grasp

# Pass the already-constructed dictionary
Grasp.run_greedy(order_dict)
