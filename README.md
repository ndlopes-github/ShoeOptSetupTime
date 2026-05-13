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
- **Heuristic Methods:** GRASP, GA and Simulated Annealing algorithms tailored to the problem structure
- **Real-World Validation:** Computational experiments based on actual production data from a footwear manufacturer demonstrate significant improvements in production efficiency

## 🚀 Solution Methods

This implementation includes:

1. **MILP (Mixed-Integer Linear Programming)** - Exact solver providing optimal solutions
2. **Split-Solve-Merge** - Novel heuristic algorithm for large-scale instances (Pg > 1)
3. **Simulated Annealing (SA)** - Tailored metaheuristic approach
4. **GRASP** - Greedy randomized adaptive search procedure tailored to the problem
5. **Genetic Algorithm (GA)** - Population-based evolutionary metaheuristic
6. **Greedy** - Simple constructive heuristic for quick solutions

## ⚡ Quick Start & Smoke Tests

### Setup

```bash
git clone https://github.com/ndlopes-github/ShoeOptSetupTime.git
cd ShoeOptSetupTime
julia --project -e 'using Pkg; Pkg.instantiate()'
```

### Verify your environment
### Inspect the scripts to check the default options

```bash
# Test 1: SSM-SA/MILP heuristic(exact) on default instance
julia --project scripts/run_ssm_milp.jl heuristic 
julia --project scripts/run_ssm_milp.jl exact 

# Test 2: Simulated Annealing on default instance
julia --project scripts/run_sa.jl 

# Test 3: GRASP on default instance
julia --project scripts/run_grasp.jl

# Test 4: Greedy on single-mold default instance
julia --project scripts/run_greedy.jl 

# Test 5: Genetic Algorithm 
julia --project scripts/run_ga.jl 
```

**Expected outcomes:**
- Solver produces solution quality metrics (cost, time, iterations)
- Output files are generated in or `data/sims/` for SSM-SA or MILP (first test), while the other output straight to terminal. 

If any test fails, ensure:
1. Julia environment is initialized: `julia --project -e 'using Pkg; Pkg.instantiate()'`
2. A solver is available: use `solver_name = "HiGHS"` in settings instance files if Gurobi is not licensed.
3. Example instance files exist in `data/settings/`. Full set of instances used in paper are available at supplementary material.

For detailed usage of each algorithm, see [Usage](#-usage).

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
Nit = 100                       # Simulated annealing iterations or Independent Runs in Grasp

                                # SA parameters
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

The `data/settings/` directory contains the real-world instances from the paper:
- **E_O2_#2_3p.jl** - Exact solver instance (Pg=1, longer time limits)
- **H_O2_#1_3p.jl** - Heuristic instance (Pg=2, single-mold, multi-stage partitioning)
- **H_O2_#2_3p.jl** - Heuristic instance (Pg=2, multi-mold, multi-stage partitioning)


### Solver Output

Each run produces:
- Log files with detailed optimization progress
- Console output showing solution quality and computation times
- CSV exports for batch experiments with statistics

## � Sample Results

Method comparison on instance **O2 / Scenario #2 / p=3** (8 jobs, 3 shelves, α=1, β=6):

| Method | Best Cost | Slots Used | Total Time (s) | Notes |
|--------|----------:|:----------:|---------------:|-------|
| MILP (Exact, Pg=1) | **1492** | 18 | 41.8 | Optimal, proven to 0% gap |
| SA (100 runs) | **1492** | 18 | 181.7 | Matches optimal |
| GA (10 runs, pop=400) | **1492** | 18 | 130.8 | Matches optimal |
| GRASP (1 run) | 1497 | 18 | 4.3 | 0.33% above optimal |
| SSM (Pg=2, Tl=30 s/sub) | 1498 | 18 | 102.8 | 0.40% above optimal |

Produced by `scripts/batch_helpers/batch_compare_all_methods.jl` on 2026-05-13 (Gurobi 13.0.0, AMD Ryzen 9 3950X).

> **SSM** uses the heuristic instance file (`H_O2_#2_3p.jl`, Pg=2, 30 s per sub-problem).
> All other methods use the exact instance file (`E_O2_#2_3p.jl`, Pg=1).

## �🔧 Troubleshooting

**Issue:** Gurobi license not found
- **Solution:** Install Gurobi or switch to HiGHS solver in settings file by changing `solver_name = "HiGHS"`

**Issue:** Out of memory errors
- **Solution:** Reduce `Tl` (time limit), decrease problem size, or increase partition size `Pg`

**Issue:** Slow performance
- **Solution:** Ensure Gurobi uses multiple threads (default 30); check system resources; try HiGHS as alternative

**Issue:** Module not found errors
- **Solution:** Ensure you are in the project directory and have run `Pkg.instantiate()` after `Pkg.activate(".")`

## 🔬 Batch Processing

For running multiple instances in batch, use the helper scripts in `scripts/batch_helpers/`. These manage sequential execution with detached screen sessions (useful for long-running experiments on remote servers).

**Available batch runners:**
- `run_ga_batch.sh` - Run GA on multiple instances
- `run_grasp_batch.sh` - Run GRASP on multiple instances  
- `run_milp_batch.sh` - Run MILP on multiple instances (WARNING: slow)
- `run_sa_batch.sh` - Run SA on multiple instances
- `run_ssm_sa_batch.sh` - Run SSM-SA on multiple instances

**Usage:**
```bash
# Run GA on instances listed in list_ga_all.txt
bash scripts/batch_helpers/run_ga_batch.sh

# Run with custom instance list
bash scripts/batch_helpers/run_grasp_batch.sh my_instances.txt

# Run with beta override
bash scripts/batch_helpers/run_sa_batch.sh scripts/batch_helpers/list_sa_all.txt --beta 6

# On remote server (survives SSH disconnect):
screen -S ga_batch -dmL -Logfile logs_ga_beta3/ga_batch.log \
    bash scripts/batch_helpers/run_ga_batch.sh
```

Results are saved to `data/exp_pro/` with logs in `logs_*_beta*/` directories.

The batch helpers are intentionally kept in `scripts/batch_helpers/` because they are maintained executable source scripts. The paper reproducibility artifacts are provided in `supplementary_material/`.

## ⚙️ Parameter Tuning with irace

For automated parameter optimization, use the irace runners in `scripts/irace_helpers/`. These are designed to work with the irace optimizer framework.

**Requirement:** irace runs require **R** with the **`irace`** package installed.

```bash
R -q -e 'install.packages("irace", repos="https://cloud.r-project.org")'
```

**Available irace runners:**
- `irace_ga_runner.jl` - GA parameter tuning target
- `irace_grasp_runner.jl` - GRASP parameter tuning target
- `irace_sa_runner.jl` - SA parameter tuning target
- `irace_ssm_sa_runner.jl` - SSM-SA parameter tuning target

**Basic usage (one parameter evaluation):**
```bash
julia --project scripts/irace_helpers/irace_sa_runner.jl \
    data/settings/H_O2_#2_3p.jl 12345 --Nit 100
```

For full irace integration, see `supplementary_material/IRACE_SCRIPTS/` for example configurations.

The `scripts/irace_helpers/` files are the maintained runners used by this codebase. The `supplementary_material/IRACE_SCRIPTS/` directory stores the paper's experiment setup and reproducibility material.

## 📊 Supplementary Material

Comprehensive experimental results, statistical analysis, and irace configuration files are available in the `supplementary_material/` directory:

- **`GA/`, `GRASP/`, `GREEDY/`, `MILP/`, `SA/`, `SSM-SA/`** - Per-method result archives (`.tar.bz2`)
- **`INSTANCES/`** - Test instances used in computational experiments
  - `INSTANCES/MILP/` - Instances for exact solver tests
  - `INSTANCES/HEURISTICS/` - Instances for heuristic method tests
- **`IRACE_SCRIPTS/`** - irace parameter tuning configurations
  - `irace_ga/`, `irace_grasp/`, `irace_sa/`, `irace_ssm_sa/` - One directory per method with scenario files, parameters, and forbidden configurations
- **`STATS_SCRIPTS/`** - Statistical analysis tools
  - `compute_rpd.jl` - Compute RPD (Relative Percent Deviation) from best known solutions
  - `statistical_analysis.jl` - Wilcoxon and Friedman tests
- **`rpd_boxplot.png`** - RPD boxplot figure
- **`supplementary_rpd.md`** - Per-instance RPD results table for all methods

## 📁 Project Structure

```
ShoeOptSetupTime/
├── scripts/                          # Algorithm implementations and runners
│   ├── run_ssm_milp.jl              # MILP and Split-Solve-Merge runner (with --help)
│   ├── split_solve_merge_milp.jl    # Novel SSM algorithm module
│   ├── run_sa.jl                    # Simulated Annealing runner (with --help)
│   ├── simulated_annealing.jl       # SA algorithm module
│   ├── run_grasp.jl                 # GRASP runner (with --help)
│   ├── grasp.jl                     # GRASP algorithm module
│   ├── run_greedy.jl                # Greedy runner
│   ├── run_ga.jl                    # Genetic Algorithm runner
│   ├── genetic_algorithm.jl         # GA algorithm module
│   ├── batch_helpers/               # Batch processing for multiple instances
│   │   ├── batch_compare_all_methods.jl  # Comprehensive comparison (with --help)
│   │   ├── run_ga_batch.sh          # GA batch sequential runner
│   │   ├── run_grasp_batch.sh       # GRASP batch sequential runner
│   │   ├── run_milp_batch.sh        # MILP batch sequential runner
│   │   ├── run_sa_batch.sh          # SA batch sequential runner
│   │   ├── run_ssm_sa_batch.sh      # SSM-SA batch sequential runner
│   │   ├── run_ga_cli.jl            # GA CLI runner
│   │   ├── run_grasp_cli.jl         # GRASP CLI runner
│   │   ├── run_milp_cli.jl          # MILP CLI runner
│   │   ├── run_sa_cli.jl            # SA CLI runner
│   │   ├── run_ssm_sa_cli.jl        # SSM-SA CLI runner
│   │   └── list_*_all.txt           # Instance lists for batch runs
│   └── irace_helpers/               # Parameter tuning runners for irace
│       ├── irace_ga_runner.jl       # GA irace target runner
│       ├── irace_grasp_runner.jl    # GRASP irace target runner
│       ├── irace_sa_runner.jl       # SA irace target runner
│       └── irace_ssm_sa_runner.jl   # SSM-SA irace target runner
├── src/
│   └── loggers.jl                   # Logging utilities
├── data/
│   ├── settings/                    # Problem instance definitions
│   │   ├── E_O2_#2_3p.jl           # Exact solver instance (Pg=1)
│   │   ├── H_O2_#1_3p.jl           # Heuristic instance (Pg=2, single-mold)
│   │   └── H_O2_#2_3p.jl           # Heuristic instance (Pg=2, multi-mold)
│   ├── exp_pro/                     # Experimental results (created at runtime)
│   └── sims/                        # Simulation logs
│       ├── exact/                   # Logs for Pg=1 runs
│       └── heuristics/              # Logs for Pg>1 runs
├── supplementary_material/          # Supplementary results and configurations
│   ├── GA/                          # GA result archives (.tar.bz2)
│   ├── GRASP/                       # GRASP result archives (.tar.bz2)
│   ├── GREEDY/                      # Greedy result archives
│   ├── MILP/                        # MILP result archives (.tar.bz2)
│   ├── SA/                          # SA result archives (.tar.bz2)
│   ├── SSM-SA/                      # SSM-SA result archives (.tar.bz2)
│   ├── INSTANCES/                   # Test instances (MILP/, HEURISTICS/)
│   ├── IRACE_SCRIPTS/               # irace tuning configurations
│   ├── STATS_SCRIPTS/               # Statistical analysis tools
│   ├── rpd_boxplot.png              # RPD boxplot figure
│   └── supplementary_rpd.md         # Per-instance RPD results table
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
