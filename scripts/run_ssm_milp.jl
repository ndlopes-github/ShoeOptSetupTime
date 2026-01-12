#= Copyright (C) 2024
Nuno David Lopes.
Created:  2025/11/21
Last changed - N. Lopes: 2026/01/12 16:20:27
Run script for Split-Solve-Merge MILP algorithm
=#

using DrWatson
@quickactivate "ShoeOptSetupTime"
using Printf

# Load settings FIRST, before including the module
include(datadir("settings", "H_O2_#2_3p.jl"))

# Now include and use the module
include(scriptsdir("split_solve_merge_milp.jl"))
using .SplitSolveMergeMILP

println("=" ^ 80)
println("Running Split-Solve-Merge MILP Algorithm")
println("=" ^ 80)

# Run SSM with timing
elapsed_time = @elapsed begin
    SplitSolveMergeMILP.run(order_dict)
end

println("\n" * "=" ^ 80)
println("SPLIT-SOLVE-MERGE EXECUTION COMPLETED")
println("=" ^ 80)
println("Total execution time: $(round(elapsed_time, digits=3))s")
println("=" ^ 80)
