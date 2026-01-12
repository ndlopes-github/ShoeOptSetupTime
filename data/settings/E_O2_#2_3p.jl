#= Copyright (C) 2024
Nuno David Lopes.
Created:  2025/11/26
Last changed - N. Lopes:2025/11/26 11:08:10
=#

using DrWatson
@quickactivate "ShoeOptSetupTime"


#######################################
# Jobs ID
g = [1 2 3 4 5 6 7 8]
# Number of molds available per job
o = [1 1 1 2 1 1 1 1]
# Job Quantities
n = [215 463 970 1240 842 342 147 99]

# Number of shelves
p = 3

# Objective function parameters
α = 1
β = 6

# Heuristics parameters
# Partition size
Pg = 1

# Number of Simulated annealing iterations
Nit = 100

# Initial temperature
T0 = 5
# Final temperature
Tf = 0.01
# Decrease (jump) iteration for temperature
Tj = 3

# Global Time Limit (seconds)
Gl = 50000

# Solver parameters
# Time limit (seconds)
Tl = 50000

# Solver selection: "Gurobi" or "HiGHS"
solver_name = "Gurobi"


#######################################
# Order Dictionary: DO NOT EDIT
order_dict = @dict g n o p α β T0 Tf Tj Pg Nit Gl Tl solver_name
# Order ID
const FILEBASENAME = splitext(basename(@__FILE__()))[1]
order_dict[:Oid] = "$(FILEBASENAME)_p_$(p)_nit_$(Nit)_Pg_$(Pg)_Tl_$(Tl)_Ts_$(T0)_$(Tf)_$(Tj)_Gl_$(Gl)"
# Add file and gitcommit stamps
@tag!(order_dict)
