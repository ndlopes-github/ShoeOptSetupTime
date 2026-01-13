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

## Quick Start: Running Individual Algorithms

The main scripts are located in the `scripts/` directory. Each algorithm has its own entry point.

### 1️⃣ MILP and Split-Solve-Merge Algorithm 

Run the Split-Solve-Merge algorithm:

```bash
# Default: runs heuristic instance (H_O2_#2_3p.jl)
julia --project scripts/run_ssm_milp.jl

# Run exact instance
julia --project scripts/run_ssm_milp.jl exact

# Run custom instance
julia --project scripts/run_ssm_milp.jl --file=E_O1_#1_2p.jl

# Show available options
julia --project scripts/run_ssm_milp.jl --help
```

The Split-Solve-Merge algorithm is a novel heuristic designed for large-scale instances. It partitions the problem, solves sub-problems independently using MILP, and merges the solutions.

### 2️⃣ Simulated Annealing

```bash
julia --project scripts/run_sa.jl
```

Runs 100 independent simulated annealing trials and reports the best solution found.

### 3️⃣ GRASP Algorithm

```bash
julia --project scripts/run_grasp.jl
```

Runs the GRASP algorithm with parameters from the settings file.

### 4️⃣ Greedy Algorithm

```bash
julia --project scripts/run_greedy.jl
```

Runs a simple greedy constructive heuristic for quick baseline solutions.

### 5️⃣ Comprehensive Batch Comparison

Compare all solution methods on all instances:

```bash
julia --project scripts/batch_compare_all_methods.jl
```

**With options:**

```bash
# Skip MILP (computationally expensive)
julia --project scripts/batch_compare_all_methods.jl --skip-milp

# Process only first 5 instances (for testing)
julia --project scripts/batch_compare_all_methods.jl --limit=5

# Skip multiple methods
julia --project scripts/batch_compare_all_methods.jl --skip-milp --skip-ssm

# Override beta parameter
julia --project scripts/batch_compare_all_methods.jl --beta=3

# Show available options
julia --project scripts/batch_compare_all_methods.jl --help
```

Output is saved to `data/exp_pro/` as CSV files with detailed comparisons.

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
- **E_O2_#2_3p.jl** - Instance used for exact/metaheuristic testing
- **H_O2_#2_3p.jl** - Instance used for heuristic algorithm testing

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
│   ├── run_milp.jl                   # Exact MILP solver
│   ├── split_solve_merge_milp.jl     # Novel Split-Solve-Merge algorithm
│   ├── run_sa.jl                     # Simulated Annealing
│   ├── run_grasp.jl                  # GRASP algorithm
│   ├── run_greedy.jl                 # Greedy heuristic
│   └── batch_compare_all_methods.jl  # Comprehensive comparison
├── src/
│   └── loggers.jl                    # Logging utilities
├── data/
│   ├── settings/                     # Problem instance definitions
│   │   ├── E_O2_#2_3p.jl            # Exact solver instance
│   │   └── H_O2_#2_3p.jl            # Heuristic solver instance
│   └── exp_pro/                      # Experimental results
├── Project.toml                      # Julia project manifest
└── README.md                         # This file
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
