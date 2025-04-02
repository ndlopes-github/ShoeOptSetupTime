#= Copyright (C) 2024
Nuno David Lopes.
Created:  2024/10/22
Last changed - N. Lopes:2025/03/25 11:20:31
=#

using DrWatson
@quickactivate "ShoeOptSetupTime"


#######################################
# EDIT 
# Jobs ID
g = [1 2 3 4 5 6 7 8]
# Number of molds available per job
o = [1 1 2 2 1 1 1 1]
# Job Quantities
n = [215 463 970 1240 842 342 147 99]

# Number of shelves 
p = 3

# Objective function parameters
α = 1
β = 3

# Heuristics parameters
# Partition size
Pg = 2

# Number of Simulated annealing iterations
Nit = 100

# Initial temperature
T0 = 5
# Final temperature
Tf = 0.01
# Decrease (jump) iteration for temperature
Tj = 3 

# Global Time Limit (seconds)
Gl= 1800

# Gurobi parameters
# Time limit (seconds)
Tl = 20
######################################

#######################################
# DO NOT EDIT #
# Order Dictionary
order_dict = @dict g n o p α β T0 Tf Tj Pg Nit Gl Tl
# Order ID
order_dict[:Oid] = "_p_$(p)_nit_$(Nit)_Pg_$(Pg)_Tl_$(Tl)_Ts_$(T0)_$(Tf)_$(Tj)_Gl_$(Gl)"
# Add file and gitcommit stamps
@tag!(order_dict)
