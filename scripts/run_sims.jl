#= Copyright (C) 2024
Nuno David Lopes.
Created:  2024/04/09
Last changed - N. Lopes: 2025/09/22 15:11:12
=#

"""
Check the simulated_annealing.jl docstrings for a detailed description.
Run the simulated annealing optimization process
run_sim(order_file="settings_file.jl")
"""

using DrWatson
@quickactivate "ShoeOptSetupTime"

include(scriptsdir("simulated_annealing.jl"))
using .SimulatedAnnealing

run_sim(; order_file="H_O2_#2_3p.jl")

run_sim(; order_file="E_O2_#2_3p.jl")

