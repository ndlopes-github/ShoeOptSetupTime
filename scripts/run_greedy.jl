#= Copyright (C) 2024
Nuno David Lopes.
Created:  2025/11/21
Last changed - N. Lopes: 2026/05/13 11:44:02
=#

using DrWatson
@quickactivate "SoftIdea"

# Load settings FIRST, before including the module
include(datadir("settings", "H_O2_#1_3p.jl"))

# Now include and use the GRASP module's greedy wrapper
include(scriptsdir("grasp.jl"))
using .Grasp

# Pass the already-constructed dictionary
Grasp.run_greedy(order_dict)
