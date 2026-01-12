#= Copyright (C) 2024
Nuno David Lopes.
Created:  2025/11/25
Run script for GRASP algorithm
=#

using DrWatson
@quickactivate "ShoeOptSetupTime"
using Printf

# Load settings
include(datadir("settings", "H_O2_#2_3p.jl"))

# Now include and use the (renamed) GRASP module
include(scriptsdir("grasp.jl"))
using .Grasp

println("=" ^ 80)
println("Running GRASP Algorithm")
println("=" ^ 80)

# Run GRASP with timing
elapsed_time = @elapsed begin
    result = Grasp.run(order_dict);
end

println("\n" * "=" ^ 80)
println("GRASP EXECUTION COMPLETED")
println("=" ^ 80)
println("Total execution time: $(round(elapsed_time, digits=3))s")
println("=" ^ 80)
