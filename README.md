# Parallel Machine Scheduling with Simultaneous Interruptions

[![Julia 1.12+](https://img.shields.io/badge/Julia-1.12+-9558B2?logo=julia&logoColor=white)](https://julialang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Status: Active](https://img.shields.io/badge/Status-Active-brightgreen.svg)](#)

This repository contains a **Julia implementation** of novel optimization algorithms for a real-world parallel machine scheduling problem in footwear manufacturing. The code accompanies a forthcoming paper and provides both exact and heuristic solution methods.

## 📄 Paper Details

**Title:** Parallel machine scheduling with simultaneous interruptions: a case study in shoe sole production

**Authors:** J. Orestes Cerdeira, Ricardo Enguiça, and Nuno Lopes

**Status:** Submitted

## 🏭 Problem Description

We address a real-world scheduling problem arising in shoe sole manufacturing, where parallel machines must be coordinated under simultaneous synchronized interruptions caused by shared resource constraints. Specifically, whenever the first job in a group of simultaneously processed jobs is completed, all remaining jobs are forced to interrupt their processing. This interruption mechanism, in which synchronized interruptions are triggered by the completion of the first job in a group (rather than by the completion of the last job or by external disruptions), leads to a scheduling structure that has not been previously studied in the literature.

> **Note:** This code is under active development. While the core algorithms are fully functional, further improvements and refinements are anticipated. We welcome contributions from the community to enhance the implementation, add features, or improve performance. Please see the [Contributing](#contributing) section for guidelines.

### 🎯 Key Contributions

- **Novel MILP Models:** Exact mixed-integer linear programming formulations of the problem
- **Split-Solve-Merge Algorithm:** A novel heuristic specifically designed for handling large-scale instances
- **Heuristic Methods:** GRASP and Simulated Annealing algorithms tailored to the problem structure
- **Real-World Validation:** Computational experiments based on actual production data from a footwear manufacturer demonstrate significant improvements in production efficiency

## 🚀 Solution Methods

This implementation includes:

1. **MILP (Mixed-Integer Linear Programming)** - Exact solver providing optimal solutions
2. **Split-Solve-Merge** - Novel heuristic algorithm for large-scale instances (Pg > 1)
3. **Simulated Annealing** - Tailored metaheuristic approach
4. **GRASP** - Greedy randomized adaptive search procedure tailored to the problem
5. **Greedy** - Simple constructive heuristic for quick solutions

## ⚡ Quick Start

### 1. Clone and Setup

```bash
git clone https://github.com/ndlopes-github/ShoeOptSetupTime.git
cd ShoeOptSetupTime
julia --project -e 'using Pkg; Pkg.instantiate()'
```

### 2. Run an Algorithm (5 minutes)

**Split-Solve-Merge (recommended for quick test):**
```bash
julia --project scripts/run_ssm_milp.jl
```

**Simulated Annealing (100 trials):**
```bash
julia --project scripts/run_sa.jl
```

**Compare All Methods:**
```bash
julia --project scripts/batch_compare_all_methods.jl --skip-milp --limit=1
```
> Note: `--skip-milp` skips the exact MILP solver as it may take a long time for larger instances. See [Usage](#-usage) section for running MILP and other options.

### 3. Custom Instance

Create `data/settings/my_instance.jl`:
```julia
using DrWatson
@quickactivate "ShoeOptSetupTime"

g = [1 2 3 4]                    # Job IDs
o = [1 1 1 1]                   # Molds available per job
n = [100 200 150 180]           # Job quantities
p = 2                           # Number of shelves
α = 1; β = 3                    # Objective coefficients
Pg = 2                          # Partition size
Tl = 20                         # Time limit (seconds)
solver_name = "Gurobi"          # or "HiGHS" for open-source

order_dict = Dict(
    :g => g, :o => o, :n => n, :p => p, :α => α, :β => β,
    :Pg => Pg, :Tl => Tl, :solver_name => solver_name
)
```

Run it:
```bash
julia --project -e 'include("data/settings/my_instance.jl"); include("scripts/split_solve_merge_milp.jl"); using .SplitSolveMergeMILP; SplitSolveMergeMILP.run(order_dict)'
```

For more detailed examples, see [Usage](#usage) section.

# 💻 Installation

## Prerequisites

- Julia 1.12 or later
- Git

## Solver Requirements

The code supports two optimization solvers:

- **Gurobi** (default, commercial license required) - Recommended for best performance
- **HiGHS** (open-source alternative) - Free and fully functional, suitable for most instances

If using HiGHS instead of Gurobi, you'll need to modify the `solver_name` parameter in the settings files from `"Gurobi"` to `"HiGHS"`.

## Setup Instructions

1. **Clone the repository:**
   ```bash
   git clone https://github.com/ndlopes-github/ShoeOptSetupTime.git
   cd ShoeOptSetupTime
   ```

2. **Set up Julia environment:**
   ```bash
   julia> using Pkg
   julia> Pkg.add("DrWatson")  # install globally for quickactivate
   julia> Pkg.activate(".")
   julia> Pkg.instantiate()
   ```

3. **Verify installation:**
   ```bash
   julia> include("Project.toml")  # loads all dependencies
   ```

This setup will install all required packages and handle local path configuration automatically.

## Dependencies

Key packages include:
- **JuMP** & **Gurobi**/**HiGHS** - Optimization modeling and solving
- **DrWatson** - Reproducible science project management
- **DataFrames** & **CSV** - Data handling and results export
- **Combinatorics** - Combinatorial operations
- **LoggingExtras** - Advanced logging capabilities

# 📖 Usage

All scripts include comprehensive built-in help and command-line options. Use `--help` to see available options for each script.

## Quick Start: Running Individual Algorithms

The main scripts are located in the `scripts/` directory. Each algorithm has its own entry point with flexible command-line options.

### 1️⃣ Exact MILP and Split-Solve-Merge Algorithm (Novel Heuristic)

The Exact MILP  and the Split-Solve-Merge algorithms are the **primary contribution** of this work.

```bash
# Show all available options
julia --project scripts/run_ssm_milp.jl --help

# Run with default heuristic instance (H_O2_#2_3p.jl)
julia --project scripts/run_ssm_milp.jl

# Or explicitly specify heuristic mode
julia --project scripts/run_ssm_milp.jl heuristic

# Run exact MILP instance (E_O2_#2_3p.jl) - takes longer
julia --project scripts/run_ssm_milp.jl exact

# Run custom instance
julia --project scripts/run_ssm_milp.jl --file=custom_instance.jl
```

**When to use:**
- Large instances where exact MILP is too slow
- When option `Pg > 1` (multi-stage partitioning)
- For production-scale problems requiring good solutions quickly

### 2️⃣ Simulated Annealing

Multiple independent SA runs provide statistical evidence of solution quality.

```bash
# Show all available options
julia --project scripts/run_sa.jl --help

# Run with defaults (100 independent runs)
julia --project scripts/run_sa.jl

# Run with fewer trials for faster testing
julia --project scripts/run_sa.jl --runs=20

# Custom instance with specific configuration
julia --project scripts/run_sa.jl --file=test_single_mold.jl --runs=50

# Silent mode (no iteration logging, only final results)
julia --project scripts/run_sa.jl --runs=10 --log-every=0

# Frequent logging (every iteration)
julia --project scripts/run_sa.jl --runs=5 --log-every=1
```

**Options:**
- `--runs=N`: Number of independent SA trials (default: 100)
- `--log-every=N`: Logging frequency per run (default: 10, use 0 for silent mode)
- `--file=PATH`: Load custom instance from `data/settings/`

### 3️⃣ GRASP Algorithm

Greedy Randomized Adaptive Search Procedure with local search.

```bash
# Show all available options
julia --project scripts/run_grasp.jl --help

# Run with default settings from instance file
julia --project scripts/run_grasp.jl

# Override iteration count from command line
julia --project scripts/run_grasp.jl --iterations=50

# Custom instance with iteration override
julia --project scripts/run_grasp.jl --file=custom_instance.jl --iterations=100
```

**Options:**
- `--iterations=N`: Override number of GRASP iterations from settings file
- `--file=PATH`: Load custom instance from `data/settings/`

### 4️⃣ Greedy Algorithm

Fast deterministic baseline heuristic - **only supports single-mold instances** (all `o[j] == 1`).

```bash
# Show all available options
julia --project scripts/run_greedy.jl --help

# Run on single-mold instance
julia --project scripts/run_greedy.jl --file=test_single_mold.jl
```

**Important:** Greedy validates input and exits with clear error if multi-mold jobs are detected.

### 5️⃣ Comprehensive Batch Comparison

Compare all solution methods across all instances with detailed CSV output.

```bash
# Show all available options
julia --project scripts/batch_compare_all_methods.jl --help

# Run all methods on all instances (WARNING: may take hours!)
julia --project scripts/batch_compare_all_methods.jl

# Skip computationally expensive exact MILP
julia --project scripts/batch_compare_all_methods.jl --skip-milp

# Process only first 5 instances (for testing)
julia --project scripts/batch_compare_all_methods.jl --limit=5

# Skip multiple methods for faster comparison
julia --project scripts/batch_compare_all_methods.jl --skip-milp --skip-ssm

# Run only SA and GRASP
julia --project scripts/batch_compare_all_methods.jl --skip-milp --skip-ssm

# Override beta parameter for all instances
julia --project scripts/batch_compare_all_methods.jl --beta=3

# Test run: 1 instance, fast methods only
julia --project scripts/batch_compare_all_methods.jl --limit=1 --skip-milp --skip-ssm
```

**Options:**
- `--limit=N, -l=N`: Process only first N instances
- `--beta=VALUE, -b=VALUE`: Override β coefficient from settings files
- `--only-file=PATH`: Run only instances listed in file (Order,Scenario,P format)
- `--skip-milp`: Skip exact MILP solver
- `--skip-ssm`: Skip Split-Solve-Merge algorithm
- `--skip-sa`: Skip Simulated Annealing
- `--skip-grasp`: Skip GRASP algorithm

**Output:** Results saved to `data/exp_pro/` as CSV files with detailed method comparisons.

## 📚 Algorithm Documentation

All algorithm modules include comprehensive docstrings following Julia conventions:

- **Function signatures** with type information
- **Algorithm descriptions** with mathematical formulations where applicable
- **Parameter explanations** with units and constraints
- **Return value specifications** with types
- **Usage examples** demonstrating typical use cases
- **Cross-references** to related functions using `@ref`

To view documentation for any function in the Julia REPL:
```julia
julia> include("scripts/split_solve_merge_milp.jl")
julia> using .SplitSolveMergeMILP
julia> ?SplitSolveMergeMILP.run
```

## ⚙️ Configuration: Creating Settings Files

Settings files define problem instances and algorithm parameters. The structure is:

```julia
# Jobs specification
g = [1 2 3 4 5 6 7 8]           # Job IDs
o = [1 1 1 2 1 1 1 1]           # Number of molds available per job
n = [215 463 970 1240 842 342 147 99]  # Job quantities

# Problem parameters
p = 3                           # Number of shelves
α = 1                           # Objective coefficient for makespan
β = 6                           # Objective coefficient for interruptions

# Algorithm-specific parameters
Pg = 2                          # Partition size for Split-Solve-Merge (Pg=1 for exact)
Nit = 100                       # Simulated annealing iterations
T0 = 5                          # Initial temperature
Tf = 0.01                       # Final temperature
Tj = 3                          # Temperature jump interval

# Solver settings
Tl = 30                         # Solver time limit (seconds)
Gl = 1800                       # Global time limit (seconds)
solver_name = "Gurobi"          # "Gurobi" or "HiGHS"
```

### Creating a Simple Custom Instance

Create a new Julia file in `data/settings/` (e.g., `custom_instance.jl`):

```julia
using DrWatson
@quickactivate "ShoeOptSetupTime"

# Simple 4-job instance
g = [1 2 3 4]                    
o = [1 1 1 1]                   
n = [100 200 150 180]            

p = 2                           
α = 1                           
β = 3                           

Pg = 1                          # Set Pg=1 for exact solution, Pg>1 for heuristic
Nit = 50                        
T0 = 2
Tf = 0.01
Tj = 2

Tl = 20                         
Gl = 300                        
solver_name = "Gurobi"          

# Create order_dict with all parameters
order_dict = Dict(
    :g => g, :o => o, :n => n, :p => p, :α => α, :β => β,
    :Pg => Pg, :Nit => Nit, :T0 => T0, :Tf => Tf, :Tj => Tj,
    :Tl => Tl, :Gl => Gl, :solver_name => solver_name
)
```

Then run any algorithm:
```bash
julia> include("data/settings/custom_instance.jl")
julia> include("scripts/split_solve_merge_milp.jl")
julia> using .SplitSolveMergeMILP
julia> SplitSolveMergeMILP.run(order_dict)
```

## 📊 Data and Instances

The `data/exp_pro/` directory contains two real-world instances from the paper:
- **E_O2_#2_3p.jl** - Exact solver instance (Pg=1, longer time limits)
- **H_O2_#2_3p.jl** - Heuristic instance (Pg=2, multi-stage partitioning)

### Creating Test Instances

For testing single-mold scenarios (required for Greedy algorithm):

```julia
# data/settings/test_single_mold.jl
using DrWatson
@quickactivate "ShoeOptSetupTime"

g = [1 2 3 4 5]
o = [1 1 1 1 1]  # All single mold
n = [100 200 150 80 120]

p = 3
α = 1
β = 6

Pg = 1
Nit = 10
Tl = 30
T0 = 5
Tf = 0.01
Tj = 3
Gl = 1800
solver_name = "Gurobi"

# Create order_dict
order_dict = @dict g n o p α β T0 Tf Tj Pg Nit Gl Tl solver_name
const FILEBASENAME = splitext(basename(@__FILE__()))[1]
order_dict[:Oid] = "$(FILEBASENAME)_p_$(p)_nit_$(Nit)_Pg_$(Pg)_Tl_$(Tl)_Ts_$(T0)_$(Tf)_$(Tj)_Gl_$(Gl)"
@tag!(order_dict)
```

Results from experiments are saved in:
- CSV comparison tables
- Detailed solver logs in separate output files
- Progress and error logs for batch runs

### Solver Output

Each run produces:
- Log files with detailed optimization progress
- Console output showing solution quality and computation times
- CSV exports for batch experiments with statistics

## 🔧 Troubleshooting

**Issue:** Gurobi license not found
- **Solution:** Install Gurobi or switch to HiGHS solver in settings file by changing `solver_name = "HiGHS"`

**Issue:** Out of memory errors
- **Solution:** Reduce `Tl` (time limit), decrease problem size, or increase partition size `Pg`

**Issue:** Slow performance
- **Solution:** Ensure Gurobi uses multiple threads (default 30); check system resources; try HiGHS as alternative

**Issue:** Module not found errors
- **Solution:** Ensure you are in the project directory and have run `Pkg.instantiate()` after `Pkg.activate(".")`

## 📁 Project Structure

```
ShoeOptSetupTime/
├── scripts/                          # Algorithm implementations
│   ├── run_ssm_milp.jl              # MILP and Split-Solve-Merge runner (with --help)
│   ├── split_solve_merge_milp.jl    # Novel SSM algorithm module
│   ├── run_sa.jl                    # Simulated Annealing runner (with --help)
│   ├── simulated_annealing.jl       # SA algorithm module
│   ├── run_grasp.jl                 # GRASP runner (with --help)
│   ├── grasp.jl                     # GRASP algorithm module
│   ├── run_greedy.jl                # Greedy runner (with --help)
│   └── batch_compare_all_methods.jl # Comprehensive comparison (with --help)
├── src/
│   └── loggers.jl                   # Logging utilities
├── data/
│   ├── settings/                    # Problem instance definitions
│   │   ├── E_O2_#2_3p.jl           # Exact solver instance (Pg=1)
│   │   ├── H_O2_#2_3p.jl           # Heuristic instance (Pg=2)
│   │   └── test_single_mold.jl     # Example single-mold instance
│   ├── exp_pro/                     # Experimental results (CSV files)
│   └── sims/                        # Simulation logs
│       ├── exact/                   # Logs for Pg=1 runs
│       └── heuristics/              # Logs for Pg>1 runs
├── Project.toml                     # Julia project dependencies
├── Manifest.toml                    # Locked dependency versions
├── LICENSE                          # MIT License
└── README.md                        # This file
```

## 📚 Citation

If you use this code in your research, please cite:

```bibtex
@article{Cerdeira2026ShoeOptSetupTime,
  title={Parallel machine scheduling with simultaneous interruptions: a case study in shoe sole production},
  author={Cerdeira, J. Orestes and Enguiça, Ricardo and Lopes, Nuno},
  journal={Submitted},
  year={2026},
  note={To appear}
}
```

## 🤝 Contributing

We welcome contributions! If you find issues or have suggestions:
1. Open an issue describing the problem
2. Submit a pull request with improvements
3. Ensure code follows existing style conventions

## 📜 License

This project is licensed under the MIT License. See the LICENSE file for details.

## 📧 Contact

For questions or inquiries about this codebase or the associated paper, please contact:
- **Nuno Lopes:** nuno(one dot)lopes(at@)isel(one dot)pt
