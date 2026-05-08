#=  Copyright (C) 2025
Nuno David Lopes.
Created: 2025/12/11
=#

module GeneticAlgorithm

export run_ga

using DrWatson
@quickactivate "SoftIdea"

using Logging, Random, Printf, Statistics, StatsBase

# Include SimulatedAnnealing module for cost function and utility functions
include(scriptsdir("simulated_annealing.jl"))
using .SimulatedAnnealing


"""
    Chromosome

Three-part chromosome encoding for the job scheduling problem with multi-mold support.

- `subjobs`: expanded list of `(jobid, moldid, quantity)` — one entry per (job, mold)
  pair. For single-mold jobs `moldid = 1`; for multi-mold jobs the quantities are
  determined at chromosome creation via `split_quantity_randomly` and can be
  redistributed by the quantity mutation operator.
- `shelf`: `shelf[i]` ∈ `1:p` — the shelf assigned to subjob `i`.
- `perm`: a permutation of `1:num_subjobs` — global priority order used to determine
  the intra-shelf ordering. Subjobs on the same shelf appear in the order they are
  encountered when iterating `perm`.
"""
mutable struct Chromosome
    subjobs::Vector{NamedTuple}
    shelf::Vector{Int}
    perm::Vector{Int}
end

"""
    expand_subjobs(g_vec, n_vec, o_vec)

Expands the job arrays into a flat list of `(jobid, moldid, quantity)` subjobs.
For single-mold jobs a single entry is produced. For multi-mold jobs,
`split_quantity_randomly` is used to split the total quantity.
"""
function expand_subjobs(g_vec::Vector, n_vec::Vector, o_vec::Vector)
    subjobs = NamedTuple[]
    for (idx, jobid) in enumerate(g_vec)
        total_qty = n_vec[idx]
        num_molds = o_vec[idx]
        if num_molds == 1
            push!(subjobs, (jobid=jobid, moldid=1, quantity=total_qty))
        else
            mold_quantities = SimulatedAnnealing.split_quantity_randomly(total_qty, num_molds)
            for (moldid, qty) in enumerate(mold_quantities)
                push!(subjobs, (jobid=jobid, moldid=moldid, quantity=qty))
            end
        end
    end
    return subjobs
end

"""
    chromosome_to_shelves(c, p)

Converts a chromosome into the `shelves` data structure used by the cost function.
Subjobs are inserted into their assigned shelf in `perm` order.
"""
function chromosome_to_shelves(c::Chromosome, p::Int)
    shelves = [Vector{NamedTuple}() for _ in 1:p]
    for subjob_idx in c.perm
        s = c.shelf[subjob_idx]
        push!(shelves[s], c.subjobs[subjob_idx])
    end
    return shelves
end

"""
    shelves_to_chromosome(shelves)

Converts a `shelves` vector back into a `Chromosome`. Subjobs are flattened in shelf
order; `perm` is set to the identity `1:n` (ordering is already encoded in the
shelf vectors). This is the inverse of `chromosome_to_shelves`.
"""
function shelves_to_chromosome(shelves::Vector)
    subjobs = NamedTuple[]
    shelf_ids = Int[]
    for (shelf_idx, shelf) in enumerate(shelves)
        for subjob in shelf
            push!(subjobs, subjob)
            push!(shelf_ids, shelf_idx)
        end
    end
    n = length(subjobs)
    return Chromosome(subjobs, shelf_ids, collect(1:n))
end

"""
    evaluate_fitness(c, p, α, β)

Calculates the fitness (cost) of a chromosome. Lower cost means better fitness.
"""
function evaluate_fitness(c::Chromosome, p::Int, α::Real, β::Real)
    shelves = chromosome_to_shelves(c, p)
    _, _, cost = SimulatedAnnealing.get_cost_from_shelves(shelves, α, β)
    return Float64(cost)
end

"""
    random_chromosome(g_vec, n_vec, o_vec, p)

Creates a random chromosome: expands jobs into subjobs (with random quantity splits
for multi-mold jobs), then assigns each subjob a random shelf and random priority.
"""
function random_chromosome(g_vec::Vector, n_vec::Vector, o_vec::Vector, p::Int)
    subjobs = expand_subjobs(g_vec, n_vec, o_vec)
    num_subjobs = length(subjobs)
    return Chromosome(subjobs, rand(1:p, num_subjobs), randperm(num_subjobs))
end

"""
    initialize_population(pop_size, g_vec, n_vec, o_vec, p)

Creates an initial population of random chromosomes.
"""
function initialize_population(pop_size::Int, g_vec::Vector, n_vec::Vector, o_vec::Vector, p::Int)
    return [random_chromosome(g_vec, n_vec, o_vec, p) for _ in 1:pop_size]
end

"""
    inverse_fitness_select(population, fitnesses, n_parents)

Selects `n_parents` individuals from `population` without replacement, with
probability inversely proportional to fitness (lower fitness = higher probability).

# Arguments
- `population`: current population vector.
- `fitnesses`: corresponding fitness values (lower is better).
- `n_parents`: number of parents to select.

# Returns
- `Vector{Chromosome}` of length `n_parents`.
"""
function inverse_fitness_select(
    population::Vector{Chromosome}, fitnesses::Vector{Float64}, n_parents::Int
)
    ε = 1e-10
    weights = Weights([1.0 / (f + ε) for f in fitnesses])
    indices = StatsBase.sample(1:length(population), weights, n_parents; replace=false)
    return population[indices]
end

"""
    ox_crossover(p1, p2)

Order Crossover (OX) for permutation vectors. Copies a random segment from `p1`
into the child, then fills the remaining positions in the order they appear in `p2`.
"""
function ox_crossover(p1::Vector{Int}, p2::Vector{Int})
    n = length(p1)
    a, b = minmax(rand(1:n), rand(1:n))
    child = fill(0, n)
    child[a:b] = p1[a:b]
    segment = Set(p1[a:b])
    remaining = [x for x in p2 if x ∉ segment]
    positions = vcat(b+1:n, 1:a-1)
    for (pos, val) in zip(positions, remaining)
        child[pos] = val
    end
    return child
end

"""
    crossover(parent1, parent2)

Performs crossover on a three-part chromosome:
- `subjobs`: child1 inherits from parent1, child2 from parent2 (quantities preserved).
- `shelf` part: single-point crossover.
- `perm` part: Order Crossover (OX), preserving permutation validity.

Returns two offspring chromosomes.
"""
function crossover(parent1::Chromosome, parent2::Chromosome)
    n = length(parent1.shelf)
    pt = rand(1:n)
    shelf1 = [parent1.shelf[1:pt]; parent2.shelf[pt+1:end]]
    shelf2 = [parent2.shelf[1:pt]; parent1.shelf[pt+1:end]]
    perm1 = ox_crossover(parent1.perm, parent2.perm)
    perm2 = ox_crossover(parent2.perm, parent1.perm)
    return Chromosome(copy(parent1.subjobs), shelf1, perm1),
    Chromosome(copy(parent2.subjobs), shelf2, perm2)
end

"""
    apply_neighbor(c, p, g, o)

Applies a SA-style `two_step_neighbor` perturbation to chromosome `c`.
Converts to shelves, applies `two_step_neighbor`, then converts back.
Returns a new `Chromosome`.
"""
function apply_neighbor(c::Chromosome, p::Int, g::Vector, o::Vector)
    shelves = chromosome_to_shelves(c, p)
    new_shelves = SimulatedAnnealing.two_step_neighbor(shelves, g, o)
    return shelves_to_chromosome(new_shelves)
end

"""
    is_clone(fitness, placed_fitnesses)

Returns `true` if `fitness` matches any value in `placed_fitnesses`
(within `rtol=1e-12`).
"""
function is_clone(fitness::Float64, placed_fitnesses::Vector{Float64})
    return any(f -> isapprox(fitness, f; rtol=1e-12), placed_fitnesses)
end

"""
    run_ga(order_dict; pop_size=100, clone_threshold=0.1, log_every=10)

Runs the genetic algorithm to find an optimal job-to-shelf assignment and intra-shelf
ordering. Supports single-mold and multi-mold jobs (`o[i] ≥ 1`).

# Arguments
- `order_dict::Dict`: Contains problem parameters `g`, `n`, `o`, `p`, `α`, `β`, `Nit`.
  `Nit` is used as the number of generations.

# Keyword Arguments
- `pop_size::Int`: Number of individuals in the population (default: 100).
- `clone_threshold::Float64`: Fraction of `pop_size` clones allowed per generation
  before excess clones are replaced by random individuals (default: 0.1).
- `log_every::Int`: Frequency of logging progress (default: 10).

# Returns
- `best_shelves`: the best partition found, as a `Vector{Vector{NamedTuple}}`.
"""
function run_ga(order_dict; pop_size::Int=100, clone_threshold::Float64=0.1,
    log_every::Int=10)
    @info "Running Genetic Algorithm with order_dict:"
    @info order_dict
    @unpack g, n, o, p, α, β, Nit = order_dict
    generations = Nit

    g_vec = vec(g)
    n_vec = vec(n)
    o_vec = vec(o)
    num_subjobs = sum(o_vec)

    max_clones = floor(Int, pop_size * clone_threshold)
    n_parents = max(1, fld(pop_size, 2))

    population = initialize_population(pop_size, g_vec, n_vec, o_vec, p)
    fitnesses = zeros(Float64, pop_size)

    best_chromosome = random_chromosome(g_vec, n_vec, o_vec, p)
    best_fitness = Inf
    best_max_sum = 0
    best_m = 0

    @info @sprintf("GA params: generations=%d (Nit) pop_size=%d n_parents=%d clone_threshold=%.3f num_subjobs=%d",
        generations, pop_size, n_parents, clone_threshold, num_subjobs)

    clone_count = 0  # clone count from the previous generation (logged each period)

    for gen in 1:generations
        # Evaluate fitness of the current population
        for i in 1:pop_size
            fitnesses[i] = evaluate_fitness(population[i], p, α, β)
        end

        # Find the best individual in the current generation
        min_fitness, min_idx = findmin(fitnesses)
        if min_fitness < best_fitness
            best_fitness = min_fitness
            best_chromosome = deepcopy(population[min_idx])
            best_shelves_temp = chromosome_to_shelves(best_chromosome, p)
            best_max_sum, best_m, _ = SimulatedAnnealing.get_cost_from_shelves(best_shelves_temp, α, β)
        end

        # Logging: gen_best = current generation's best; best_cost = running overall best
        if log_every > 0 && gen % log_every == 0
            @info @sprintf("gen=%5d gen_best=%.6f best_cost=%.6f best_m=%d best_max_sum=%.6f clones=%d",
                gen, min_fitness, best_fitness, best_m, best_max_sum, clone_count)
        end

        # Select n_parents via inverse-fitness weighting
        parents = inverse_fitness_select(population, fitnesses, n_parents)

        # Generate 2*n_parents offspring; clone check against current population + placed children
        placed_fitnesses = copy(fitnesses)  # includes current population (paper §2.4 clone scope)
        offspring = Chromosome[]
        offspring_fitnesses = Float64[]
        clone_count = 0

        for _ in 1:n_parents
            pa = parents[rand(1:n_parents)]
            pb = parents[rand(1:n_parents)]
            child1, child2 = crossover(pa, pb)
            child1 = apply_neighbor(child1, p, g_vec, o_vec)
            child2 = apply_neighbor(child2, p, g_vec, o_vec)
            for child in (child1, child2)
                f = evaluate_fitness(child, p, α, β)
                if is_clone(f, placed_fitnesses) && clone_count >= max_clones
                    child = random_chromosome(g_vec, n_vec, o_vec, p)
                    f = evaluate_fitness(child, p, α, β)
                else
                    is_clone(f, placed_fitnesses) && (clone_count += 1)
                end
                push!(offspring, child)
                push!(offspring_fitnesses, f)
                push!(placed_fitnesses, f)
            end
        end

        # Merge old population + offspring; keep best pop_size (paper §2.4)
        combined = [population; offspring]
        combined_fit = [fitnesses; offspring_fitnesses]
        best_indices = sortperm(combined_fit)[1:pop_size]
        population = combined[best_indices]
    end

    @info @sprintf("GA finished: generations=%d best_cost=%.6f best_m=%d best_max_sum=%.6f",
        generations, best_fitness, best_m, best_max_sum)

    best_shelves = chromosome_to_shelves(best_chromosome, p)

    @info "Best partition (shelves) snapshot:"
    for (i, shelf) in enumerate(best_shelves)
        println("  Shelf $i: ", join(["(job=$(j.jobid), mold=$(j.moldid), qty=$(j.quantity))" for j in shelf], ", "))
    end

    df = SimulatedAnnealing.partition_to_dataframe(best_shelves)
    @info df
    dfp = SimulatedAnnealing.extract_slots_from_partition(df)
    @info dfp

    return best_shelves
end

end # module GeneticAlgorithm
