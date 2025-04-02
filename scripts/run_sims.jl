#= Copyright (C) 2024
Nuno David Lopes.
Created:  2024/04/09
Last changed - N. Lopes:2025/04/01 15:37:27
=#
using DrWatson
@quickactivate "ShoeOptSetupTime"

include(scriptsdir("simulated_annealing.jl"))
using .SimulatedAnnealing

"""
Check the simulated_annealing.jl docstrings for a detailed description.
Run the simulated annealing optimization process
run_sim(order_file="settings_file.jl")
"""
run_sim(; order_file="H_O1_33.jl")
