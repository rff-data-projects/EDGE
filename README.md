# EDGE
The Employment Dynamics in General Equilibrium model is an extension of modeling work published in Hafstead and Williams (2018) and Hafstead et al. (2022). EDGE is used to evaluate how environmental policies such as carbon pricing impact unemployment and the reallocation of workers across sectors.

A complete description of the benchmark model and its extensions is available [here](https://media.rff.org/documents/EDGE_Model_Documentation-070122.pdf).

EDGE uses GAMS with the PATH solver.

Each folder contains a variation of the EDGE model and new variations will be added over time as they become available.

Within each model folder, the "MAIN" files are completely self-contained files to simulate the effect of carbon pricing on the model economy  They define primary or "deep" parameters, solves for a benchmark steady state, defines the model, defines price paths, verifies the calibration and policy models are consistent, and runs a pre-specified number of policy simulations.
