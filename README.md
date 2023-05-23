# Ising1d

Simulate the 1-dimensional classical Ising model with periodic boundary
conditions using Gibbs sampling.

## Basic Usage

Start a simulation by running the command:

    ./ising1d path/to/input_file.toml

where input_file.toml is a TOML file containing model parameters. An example
input file might look like:

    sites = 1024
    init_iters = 10000
    sample_iters = 100000
    sample_freq = 100
    coupling = 1.5
    field = 0.0

which would correspond to a simulation of 1024 ising spins with an interaction
between sites of J/kBT = 1.5, an external electric field strength of h/kBT = 0.0. 
In this example 10,000 sample configurations are drawn using Gibbs sampling and
discarded to allow the system to equilibrate and then 100,000 samples are
subsequently drawn. Of these 100,000 sample configurations, every 100 of them
are used to estimate the net magnetization of the ensemble as well as
inter-site correlations.
