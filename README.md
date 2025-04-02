# ShoeOptSetupTime 

This repository implements an optimization algorithm based on a forthcoming paper. The contents of this repository are under active development and may undergo changes to align with the paper's final version. Please stay tuned for updates as we continue refining the implementation.

## Paper Details

   > Title: 

   > Authors: [To appear]

<!---
   > Authors: J. O. Cerdeira 1, R. Enguiça 2,  N. Lopes 3, A. Moura 4

   > 1- CMA, Department of Mathematics, NOVA University Lisbon; 
     2- ISEL, Polytechnic of Lisboa, and CEMAT, University of Lisboa;
     3- ISEL, Polytechnic of Lisboa, and CEMAT, University of Lisboa;
     4- ISEP-LEMA, Polytechnic of Porto, and CMUP, University of Porto;
--->
   > Journal: [To appear]

   > Publication Date: [To appear]

   > Abstract: [To appear]


# Description

The code in this repository aims to replicate the algorithm described in the forthcoming paper. Our intention is to provide a practical implementation that can be used and tested by the community. As such, please consider this code as a work in progress, subject to further modifications and improvements.


This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> ShoeOptSetupTime

It is authored by N.Lopes.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "ShoeOptSetupTime"
```
which auto-activate the project and enable local path handling from DrWatson.

# USAGE


## Running the Model

To run the model, follow the steps below:

1. Open a console or terminal at the project directory.

2. Run the following command:
   ```
   $julia scripts/runsims.jl
   ```
      >   This command executes an instance described in the settings file "H_O1_33.jl" in data/settings dir.

## Settings File

The parameters are the following
```
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
```

### Output

After running the model, you can expect the following output:

+ An info log.txt file reflecting the heuristic or optimal solution for the instance.

### Troubleshooting

    Ensure the settings file is correctly formatted and all required parameters are provided.

    Verify that the required Julia packages (DrWatson, Gurobi, etc.) are installed and accessible.

# Contributing

We welcome contributions to this project! If you find any issues or have suggestions for improvements, please feel free to open an issue or submit a pull request. We appreciate your involvement in making this implementation more robust and accurate.

# License

This project is currently under MIT license. Please refer to the LICENSE file for more information.

# Contact

If you have any questions or inquiries regarding this codebase or the associated paper, please contact: [nuno(dot)lopes(at)isel(dot)pt]
