# Lattice Gas Simulation in C

This project implements a 2D lattice gas model in C, inspired by the Ising model and using Monte Carlo techniques. The goal is to study phase transitions in systems of interacting particles.

## 🧠 Description

Each site on a 2D square lattice is either occupied or empty (binary state). The system evolves under fixed temperature and particle number (canonical ensemble), using an energy-based probabilistic update rule. The interaction energy depends on neighboring particles, and boundary conditions are periodic.

The model is mathematically related to the 2D Ising model via a variable transformation.

## 🔧 Features

- Monte Carlo simulation of particle hopping
- Periodic boundary conditions
- Energy calculation using local neighbor interactions
- Temperature-dependent phase behavior (gas / condensed phases)
- Automatic generation of energy curves and grid visualizations via Gnuplot
- Determination of critical temperature `Tc` as a function of particle density
- Study of:
  - Temporal correlation functions
  - Spatial correlation functions
  - Cluster size (coherence length)
  - Influence of negative interaction constants (`J < 0`)

## Requirements

- **C compiler** (e.g., `gcc`)
- **Gnuplot** for visualization

