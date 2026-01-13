# Test instance with all single-mold jobs for greedy algorithm testing

using DrWatson
@quickactivate "ShoeOptSetupTime"

# Jobs ID
g = [1 2 3 4 5]
# Number of molds available per job (all single mold)
o = [1 1 1 1 1]
# Job Quantities
n = [100 200 150 80 120]

# Number of shelves
p = 3

# Objective function parameters
α = 1
β = 6

# Heuristics parameters
Pg = 1
Nit = 10

# Solver time limit (seconds)
Tl = 30

# Solver Name
solver_name = "Gurobi"

# Simulated annealing parameters
T0 = 5      # Initial temperature
Tf = 0.01   # Final temperature
Tj = 3      # Jump parameter
Gl = 1800   # Global time limit for SA

# Order Dictionary: DO NOT EDIT
order_dict = @dict g n o p α β T0 Tf Tj Pg Nit Gl Tl solver_name
# Order ID
const FILEBASENAME = splitext(basename(@__FILE__()))[1]
order_dict[:Oid] = "$(FILEBASENAME)_p_$(p)_nit_$(Nit)_Pg_$(Pg)_Tl_$(Tl)_Ts_$(T0)_$(Tf)_$(Tj)_Gl_$(Gl)"
# Add file and gitcommit stamps
@tag!(order_dict)
