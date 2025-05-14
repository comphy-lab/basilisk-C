# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Structure

This repository contains Basilisk, a high-performance computational fluid dynamics (CFD) framework:
- `basilisk-source/`: Official source code with the core solvers and functionality
- `basilisk-wiki/`: Documentation and sandbox implementations where researchers experiment

## Basilisk Scientific Framework

### Core Numerical Methods

- **Finite Volume Methods**: Implements conservative schemes on adaptive Cartesian grids
- **Volume-of-Fluid (VOF)**: Used for interface tracking in multiphase flows (`vof.h`)
- **Adaptive Mesh Refinement**: Quad/octree-based mesh adaptation for optimal resolution
- **Time Integration**: Various schemes from explicit to implicit methods (RK, predictor-corrector)
- **Projection Methods**: For enforcing incompressibility in Navier-Stokes solutions

### Physical Models

- **Navier-Stokes Equations**: Core incompressible fluid solver (`navier-stokes/*.h`)
- **Multiphase Flows**: Two-phase implementations with surface tension (`two-phase.h`)
- **Free Surface Flows**: Saint-Venant equations for shallow water (`saint-venant.h`)
- **Compressible Flows**: Solvers for compressible fluid dynamics (`compressible.h`)
- **Surface Tension**: Accurate methods for interface tension and curvature calculation

### Coordinate Systems

- **Cartesian**: Default coordinate system
- **Axisymmetric**: For problems with rotational symmetry (`axi.h`)
- **Spherical**: For geophysical and astrophysical applications (`spherical.h`)
- **General Orthogonal**: For complex domain representation

## Code Knowledge

### Basilisk C Language Extensions

Basilisk extends C with domain-specific language elements for PDEs:
- **Scalar/Vector Fields**: First-class data types for physical quantities
- **Foreach Loops**: Iteration over grid cells with built-in parallelism
- **Stencil Operators**: High-level gradient, interpolation operators
- **Event System**: For scheduling code execution at specific times
- **Attributes**: For defining field properties and boundary conditions

### Key Headers and Their Functions

- `grid/*`: Grid implementation (Cartesian, tree, octree) with data structures
- `common.h`: Common definitions, initialization, and finalization
- `run.h`: Setup for simulation runs, timestepping, and events system
- `poisson.h`: Multigrid Poisson solver (pressure projection)
- `diffusion.h`: Solvers for diffusion terms
- `vof.h`: Volume-of-Fluid interface tracking
- `two-phase.h`: Two-phase flow implementation
- `tension.h`: Surface tension models
- `curvature.h`: Interface curvature calculation

### Visualization System

- `output.h`: File output in various formats
- `view.h`: Interactive visualization capabilities
- `bview.c`: The 2D/3D viewer implementation
- Functions like `output_ppm()`, `output_gfs()`, and `output_vtk()`

## Analysis Capabilities

When examining Basilisk simulations, Claude can help with:

1. **Stability Analysis**: Identifying potential CFL violations, numerical instabilities
2. **Conservation Properties**: Verifying mass/momentum conservation
3. **Numerical Accuracy**: Assessing discretization schemes and convergence
4. **Boundary Conditions**: Implementing and analyzing boundary treatments
5. **Parameter Studies**: Systematically varying parameters for sensitivity analysis

## Example Assistance Tasks

Claude can help with tasks like:

- Explaining complex parts of solvers (e.g., pressure projection methods)
- Identifying appropriate numerical schemes for specific physical problems
- Developing appropriate boundary conditions for complex geometries
- Implementing new physical models on top of existing framework
- Diagnosing convergence or stability issues in simulations
- Visualizing and analyzing simulation results

## Command Reference

Essential commands for working with the codebase:
- Compile: `qcc -o program program.c`
- Run with visualization: `./program && bview2D -s program.plt`
- Run tests: `./runtest test.c`