# Molecular-Dynamics-Simulator

A C++ code for molecular dynamics simulations using a Lennard-Jones potential and a Velocity-Verlet integrator. In the future, more potentials and integrators will be added for benchmarking and comparisons.

## Prerequisites

- A C++17 compiler (e.g. `g++`, `clang++`).
- CMake (version ≥ 3.10).
- OpenMP.
- Gnuplot for data visualization.
- VMD for visualizing `.pdb` files.

## Build

```bash
mkdir build && cd build
cmake ..
make
```

## Execution

```bash
./md_simulator
```

This will run the simulation, produce benchmark outputs, and generate a `.pdb` file.

## Visualizing with VMD

You can visualize `out.pdb` in VMD:

```bash
vmd out.pdb
```

Then go to `Graphics`, press `Representation` then select `Licorice`

## Data

Simulation metrics are stored in `.dat` files under the `data/` folder. These files can be plotted with Gnuplot or any other plotting tool.

## Gnuplot

Several Gnuplot scripts are provided in the root folder of the project to plot the data found in `data/`.

## Benchmark

The code includes a `Benchmark` class for measuring iteration times and counting approximate operations (ops). **Note**: The GFLOPS count is only approximate and may be incorrect, as it is based on manual instrumentation in each function. Don't use it neither for relative comparisons or absolute metric since I didn't had time to finish my wrapper (A C++ reader to count automatically Gflop/s by detecting operation etc.)

## Author

**Msilini Yassine**  

