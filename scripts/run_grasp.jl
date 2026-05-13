#= Copyright (C) 2024
Nuno David Lopes.
Created:  2025/11/25
Run script for GRASP algorithm
=#

using DrWatson
@quickactivate "SoftIdea"
using Printf

# Beta override: set to nothing to use the value from the instance file,
# or set a number (e.g. 12) to bypass it.
const BETA_OVERRIDE = nothing

# Load settings
include(datadir("settings", "H_O2_#2_3p.jl"))

if !isnothing(BETA_OVERRIDE)
    order_dict[:β] = BETA_OVERRIDE
    println("Beta override active: β = $BETA_OVERRIDE (instance value bypassed)")
end

# Now include and use the (renamed) GRASP module
include(scriptsdir("grasp.jl"))
using .Grasp

println("="^80)
println("Running GRASP Algorithm")
println("="^80)

# Run GRASP with timing
# HACK TO REMOVE OVERHEAD COMPILATION TIMES
Grasp.run(order_dict);

elapsed_time = @elapsed begin
    result = Grasp.run(order_dict)
end

println("\n" * "="^80)
println("GRASP EXECUTION COMPLETED")
println("="^80)
println("Total execution time: $(round(elapsed_time, digits=3))s")
println("="^80)
