# Ising1d

Simulate the 1-dimensional classical Ising model with periodic boundary
conditions using Gibbs sampling.

## Basic Usage

Start a simulation by running the command:

    ./ising1d path/to/input_file.toml

where input_file.toml is a TOML file containing model parameters. An example
input file might look like:

    sites = 1024                    # Number of spin sites
    init_iters = 10000              # Number of Gibbs samples to draw to burn in the simulation
    sample_iters = 100000           # Number of Gibbs samples to draw during the main simulation
    sample_freq = 100               # During production, how frequently to collect observables
    coupling = 1.5                  # Spin-spin interaction strength, J / kB*T
    field = 0.0                     # Magnetic field strength, h / kB*T
    correlation_file = "corr.dat"   # Output file for inter-site correlations
    seed = 42                       # Random seed, for reproducibility. Omit to generate a seed randomly
