# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Structure

This repository contains the comphy-lab fork of Basilisk, a high-performance computational fluid dynamics (CFD) framework:
- `basilisk-source/`: Official source code with the core solvers and functionality (synced from Darcs)
- `basilisk-wiki/`: Documentation and wiki content (synced from Darcs)
- `bugs/`: Test cases and minimal reproductions for tracking bugs and issues
- Both main directories are synchronized daily via GitHub Actions with http://basilisk.fr

## Build System and Development Commands

### Installation and Setup
```bash
# Install dependencies (Debian/Ubuntu)
sudo apt install darcs make gawk gnuplot imagemagick ffmpeg graphviz valgrind

# Set up Basilisk
cd basilisk-source/src
ln -s config.gcc config    # or config.clang, config.osx, etc.
make
export BASILISK=$PWD
export PATH=$PATH:$BASILISK
```

### Core Development Commands
- **Compile program**: `qcc -o program program.c` (Basilisk's custom compiler)
- **Compile with options**: `qcc -O2 -Wall -disable-dimensions program.c -o program -lm`
- **Run test**: `./runtest test.c` (automated test runner with 3-hour timeout)
- **Clean build**: `make -k clean && make`
- **Update from upstream**: `darcs pull && make -k clean && make`

### Testing and Debugging
- **Single test**: `qcc test.c -o test && ./test`
- **Test with visualization**: `./program && bview2D -s program.plt`
- **Memory debugging**: Use built-in memory tracing (see README.mtrace)
- **Parallel debugging**: Built-in MPI support and profiling tools

### Configuration Files
Choose appropriate config file based on your system:
- `config.gcc`: Standard GCC configuration
- `config.clang`: Clang compiler
- `config.osx`: macOS-specific settings
- `config.gcc.32bits`: 32-bit systems

## Basilisk Scientific Framework Architecture

### High-Level Design Philosophy
Basilisk extends C with domain-specific language elements for solving PDEs on adaptive Cartesian meshes. The framework uses a literate programming approach where documentation and code coexist.

### Core System Components

**Event System** (`grid/events.h`):
- Schedules code execution at specific times/iterations
- Automatically classifies expressions as initialization, condition, or increment
- **Known Issue**: Floating-point precision bug for time increments < 1.2e-7 (see `bugs/` directory)

**Grid Management** (`grid/`):
- Adaptive mesh refinement using quad/octree structures
- Support for Cartesian, axisymmetric, and spherical coordinates
- Built-in parallelization through foreach loops

**Solver Architecture**:
- **Projection Methods**: Navier-Stokes solver with pressure projection (`poisson.h`)
- **Interface Tracking**: Volume-of-Fluid (VOF) method (`vof.h`)
- **Multiphase Flows**: Two-phase implementations with surface tension (`two-phase.h`)
- **Time Integration**: Multiple schemes from explicit to implicit (`predictor-corrector.h`, `runge-kutta.h`)

### Language Extensions
Basilisk adds these first-class constructs to C:
- **Scalar/Vector Fields**: Physical quantities with automatic memory management
- **Foreach Loops**: Grid iteration with built-in parallelism
- **Stencil Operators**: High-level gradient and interpolation operators
- **Events**: Time/iteration-based code execution scheduling
- **Attributes**: Field properties and boundary condition definitions

### Visualization System
- **Interactive 3D**: WebGL-based visualization through `bview.c`
- **Output Formats**: PPM, GFS, VTK support through `output.h`
- **Animation Tools**: Built-in support for creating videos (ppm2mp4, ppm2mpeg)

## Development Workflow

### Bug Tracking and Testing
This fork specifically tracks comphy-lab related issues. When creating test cases:
1. Place minimal reproductions in `bugs/` directory
2. Use naming convention: `issue-description.c` and `runTest_issue-description.sh`
3. Include both working and failing cases for comparison
4. Document expected vs actual behavior

### Integration with Main Basilisk
- This repository syncs daily with basilisk.fr Darcs repositories
- Maintain compatibility with upstream changes
- Report bugs both here (GitHub Issues) and to [Basilisk Google Group](https://groups.google.com/d/forum/basilisk-fr)

### Code Analysis Capabilities
When examining Basilisk simulations, focus on:
- **Numerical Stability**: CFL conditions, convergence analysis
- **Conservation Properties**: Mass, momentum, energy conservation
- **Boundary Conditions**: Proper implementation for complex geometries
- **Memory Management**: Basilisk's automatic grid adaptation and field management
- **Performance**: Parallel efficiency and memory usage

## Key Headers and Dependencies

**Essential Headers**:
- `common.h`: Core definitions and initialization
- `run.h`: Simulation setup and event system
- `navier-stokes/centered.h`: Main flow solver
- `grid/multigrid.h`: Grid and multigrid solver
- `view.h`: Visualization capabilities

**Solver-Specific Headers**:
- `saint-venant.h`: Shallow water equations
- `compressible.h`: Compressible flow dynamics
- `tension.h`, `curvature.h`: Surface tension modeling
- `diffusion.h`: Diffusion solvers
- `embed.h`: Embedded boundary methods

## Installation Notes

### System Requirements
- C99-compliant compiler (GCC recommended)
- GNU make or compatible
- Darcs for version control (install via Cabal if not in package repos)

### Optional Dependencies
- **Python Integration**: Requires SWIG and python-dev
- **Advanced Visualization**: gnuplot, imagemagick, ffmpeg
- **Fortran Components**: gfortran for GOTM, CVMix, PPR modules

### Platform-Specific Notes
- **macOS**: Use `config.osx` and see basilisk-wiki for macOS installation guide
- **32-bit Systems**: Use `config.gcc.32bits`
- **Ubuntu 24.04**: Special installation notes for OpenConnect VPN in Tips.md