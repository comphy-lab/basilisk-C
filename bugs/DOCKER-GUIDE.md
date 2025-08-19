# Basilisk Docker Development Guide

Complete guide for developing Basilisk simulations on macOS using Docker. This unified guide provides everything you need to compile, run, and debug master-slave coupling simulations.

## 🚀 Quick Start (5 Minutes)

**First time setup:**
```bash
# 1. Build and start Docker environment
./run_docker.sh

# 2. Access container and compile
./docker_shell.sh
cd /host
./compile.sh

# 3. Run simulation
./master-mac

# 4. Copy results (exit container first)
exit
./copy_from_docker.sh
```

**Daily workflow:**
```bash
# Edit files on Mac (any editor)
vim master-mac.c

# Compile in Docker
./docker_shell.sh
cd /host && ./compile.sh && ./master-mac
exit

# Get results
./copy_from_docker.sh
```

## 🎯 Overview

### Why Docker?
The Basilisk master-slave coupling requires GNU `objcopy` with symbol filtering capabilities that are **not available on macOS**. This Docker setup provides a complete Ubuntu Linux environment with all necessary GNU toolchain components.

### What This Solves
- ❌ macOS lacks GNU `objcopy` for symbol filtering
- ❌ Symbol conflicts prevent master-slave compilation
- ❌ BSD toolchain incompatibilities with Basilisk build system
- ✅ Complete GNU Linux environment in Docker
- ✅ Seamless file sharing between Mac and container
- ✅ Full Basilisk compatibility with all features

## 📋 Prerequisites

- **Docker Desktop** installed and running on macOS
- **basilisk-C repository** cloned locally
- **Terminal** with bash/zsh shell

## 🛠️ Installation & Setup

### Initial Setup
1. **Install Docker Desktop** from [docker.com](https://docker.com)
2. **Start Docker** (whale icon in menu bar should show "Docker Desktop is running")
3. **Navigate to bugs directory:**
   ```bash
   cd basilisk-C/bugs
   ```

### Build Docker Environment
```bash
./run_docker.sh
```

This builds an Ubuntu 22.04 container with:
- All Basilisk dependencies (gcc, make, gawk, gnuplot, etc.)
- GNU binutils (including `objcopy`)
- Development tools (vim, git, valgrind)
- Your local directory mounted as `/host`

## 🔄 Complete Development Workflow

### Core Principle: Edit Locally, Build in Docker

**File Organization:**
- 🍎 **Mac**: Edit source code (`master-mac.c`, `slave-mac.c`)
- 🐳 **Docker**: Compile and run (`qcc`, `objcopy`, execution)
- 📁 **Results**: Available both locally and via copy script

### Step-by-Step Process

#### 1. Edit Source Code (Mac)
```bash
# Use any Mac editor
vim master-mac.c
code slave-mac.c
nano master-mac.c
```

#### 2. Compile in Docker
```bash
# Enter Docker container
./docker_shell.sh

# Navigate to mounted directory
cd /host                # This is your Mac's bugs/ directory

# Compile the simulation
./compile.sh            # Complete build script
```

#### 3. Run Simulation (Docker)
```bash
# Still in /host directory
./master-mac            # Execute coupled simulation
```

#### 4. Access Results (Mac)
```bash
# Exit Docker container
exit

# Option A: Results already visible locally
ls -la *.mp4 master-mac

# Option B: Organized copying
./copy_from_docker.sh   # Copies to docker-results/
```

### Iteration Cycle
```bash
# Fast iteration (container stays open)
./docker_shell.sh
cd /host

# Repeat as needed:
./compile.sh && ./master-mac
./compile.sh && timeout 30s ./master-mac
```

## 📁 File Organization & Purposes

### Essential Files
| File | Purpose | Location | Usage |
|------|---------|----------|--------|
| `Dockerfile` | Ubuntu container definition | Mac | Auto-used by run_docker.sh |
| `run_docker.sh` | Build & launch container | Mac | First-time setup |
| `docker_shell.sh` | Multi-purpose container access | Mac | Daily development |
| `compile.sh` | Complete compilation script | Both | Build simulations |
| `copy_from_docker.sh` | Extract results to local | Mac | Get organized output |

### Source Code
| File | Purpose | Edit Location |
|------|---------|---------------|
| `master-mac.c` | Master simulation (channel flow) | Mac (any editor) |
| `slave-mac.c` | Slave simulation (cylinder flow) | Mac (any editor) |

### Generated Files
| File | Created By | Purpose |
|------|------------|---------|
| `master-mac` | compile.sh | Executable simulation |
| `slave.o` | compile.sh | Filtered object file |
| `vort-master.mp4` | master-mac | Master vorticity animation |
| `vort-slave.mp4` | master-mac | Slave vorticity animation |

### Directory Structure
```
basilisk-C/bugs/                    # 🍎 Mac (main workspace)
├── master-mac.c                    # ✏️  Edit here
├── slave-mac.c                     # ✏️  Edit here
├── compile.sh                      # 🔨 Build script
├── docker_shell.sh                # 🐳 Docker access
├── run_docker.sh                  # 🐳 Docker setup
├── copy_from_docker.sh            # 📁 File management
├── Dockerfile                      # 🐳 Container definition
├── docker-results/                # 📊 Organized results
│   ├── master-mac                  # Executable
│   ├── *.mp4                      # Movies
│   └── *.dat                      # Data files
└── [build outputs]                # Direct results

Container View:
├── /host/                          # 🔗 Your Mac files (use this!)
├── /workspace/basilisk-C/          # 🐳 Container copy (don't edit)
└── /workspace/basilisk-C/bugs/     # 🐳 Container copy (don't edit)
```

### What NOT to Do
- ❌ Don't edit files in `/workspace/basilisk-C/bugs/` (container copy)
- ❌ Don't manually sync files between locations
- ❌ Don't try to compile on Mac (missing objcopy)

### What TO Do
- ✅ Edit all files on Mac in your `bugs/` directory
- ✅ Use `/host/` directory when inside Docker container
- ✅ Let Docker mounting handle file sharing automatically

## 🔧 Scripts Reference

### `run_docker.sh` - Initial Setup
**Purpose:** Builds and launches the Docker container.

```bash
./run_docker.sh
```

**What it does:**
- Builds Docker image with Ubuntu 22.04 + all dependencies
- Mounts current directory as `/host`
- Mounts parent directory as `/project`
- Sets up Basilisk environment variables
- Provides interactive shell access

---

### `docker_shell.sh` - Multi-Purpose Container Access
**Purpose:** Flexible container interaction script.

**Usage:**
```bash
./docker_shell.sh [command]
```

**Commands:**
- `shell` or no argument: Open interactive bash shell
- `run`: Execute master-mac simulation
- `compile`: Run compilation script
- `status`: Show container and file status
- `stop`: Stop the container
- `start`: Start stopped container
- `logs`: Show container logs

**Examples:**
```bash
./docker_shell.sh              # Most common - open shell
./docker_shell.sh run          # Quick simulation run
./docker_shell.sh compile      # Quick recompile
./docker_shell.sh status       # Check what's running
```

---

### `compile.sh` - Complete Compilation Script
**Location:** Works in both Mac and Docker `/host/` directory

**Purpose:** Handles complete master-slave compilation with proper symbol filtering.

```bash
# Inside Docker container
cd /host
./compile.sh
```

**What it does:**
1. **Environment Setup**: Sets `BASILISK` path and environment variables
2. **File Validation**: Checks for required source files (`master-mac.c`, `slave-mac.c`)
3. **Clean Build**: Removes previous build artifacts
4. **Slave Compilation**: 
   ```bash
   qcc -D_OBJECT -c slave-mac.c -o slave.o
   ```
5. **Symbol Filtering**: Uses GNU `objcopy` to keep only coupling functions:
   - `slave_step`
   - `slave_stop`
   - `slave_interpolate`
6. **Master Linking**: Links master with filtered slave object
7. **Verification**: Confirms successful build

**Output Files:**
- `slave.o` - Filtered object file (Linux format)
- `master-mac` - Executable coupled simulation

---

### `copy_from_docker.sh` - Result Management
**Purpose:** Organized copying of simulation results from container to Mac.

```bash
./copy_from_docker.sh [type]
```

**Copy Types:**
- `all` (default): Copy all simulation results
- `executable` or `exe`: Copy just the master-mac executable
- `movies`: Copy animation files (`*.mp4`, `*.ppm`)
- `data`: Copy data files (`*.dat`, `*.log`, `*.txt`)
- `list`: List available files in container

**Examples:**
```bash
./copy_from_docker.sh           # Copy everything
./copy_from_docker.sh exe       # Just executable
./copy_from_docker.sh movies    # Just animations
./copy_from_docker.sh list      # See what's available
```

**Output:** All files copied to `./docker-results/` directory.

## 🐳 Docker Environment Details

### Container Specifications
- **Base Image:** Ubuntu 22.04 LTS
- **Architecture:** Multi-arch (ARM64 for Apple Silicon, AMD64 for Intel)
- **Container Name:** `basilisk-work`
- **Auto-restart:** Yes (unless explicitly stopped)

### Installed Packages
**Essential Basilisk Dependencies:**
```
build-essential     # gcc, make, g++, etc.
gawk               # GNU awk (required by Basilisk)
gnuplot            # Plotting and visualization
imagemagick        # Image processing
ffmpeg             # Video generation
graphviz           # Graph visualization
valgrind           # Memory debugging
git                # Version control
darcs              # Basilisk's version control
binutils           # GNU objcopy (key for symbol filtering)
vim                # Text editor
```

### Directory Mapping
| Host Location | Container Location | Purpose |
|---------------|-------------------|---------|
| `bugs/` | `/host` | **Primary workspace** (read/write) |
| `basilisk-C/` | `/project` | Parent directory access |
| (container) | `/workspace/basilisk-C/` | Basilisk repository copy |

### Environment Variables (Container)
```bash
BASILISK=/workspace/basilisk-C/basilisk-source/src
PATH=$PATH:$BASILISK
```

## 🎯 Simulation Details

### Master-Slave Coupling Architecture

The simulation demonstrates Basilisk's capability to couple two different fluid solvers in real-time:

#### Slave Simulation (Cylinder Flow)
```c
// slave-mac.c
```
- **Solver:** Navier-Stokes with embedded boundaries
- **Geometry:** Flow around circular cylinder
- **Physics:** von Kármán vortex street formation
- **Domain:** 8×8 physical units, 128×128 grid resolution
- **Time Phase:** t = 0 → 15 seconds (initialization)
- **Output:** Velocity field data for coupling

#### Master Simulation (Channel Flow)
```c
// master-mac.c
```
- **Solver:** Navier-Stokes without embedded boundaries
- **Geometry:** Straight channel flow
- **Boundary Condition:** Left inlet uses slave velocity field
- **Time Phase:** t = 15 → 30 seconds (coupled operation)
- **Output:** Combined flow visualization

#### Coupling Mechanism
- **Function:** `slave_interpolate()` - Samples slave velocity field
- **Location:** Along vertical line at slave domain exit
- **Frequency:** Every time step during master simulation
- **Data Transfer:** Real-time velocity components (u, v)
- **Synchronization:** Master time stepping controls slave queries

### Expected Simulation Behavior

#### Console Output Pattern
```
slave 0 0 0.00454545 0 0
slave 1 0.00454545 0.00454545 3 1
slave 2 0.00909091 0.00454545 5 1
...
slave 3300 15 0.00454545 2 1
master 0 0 0.0039416 0 0
master 1 0.0039416 0.0039416 2 1
master 2 0.007883 0.0039416 4 1
...
```

**Output Format:** `simulation step time dt mgp mgi`
- `step`: Time step number
- `time`: Current simulation time
- `dt`: Time step size
- `mgp`: Multigrid P-cycles
- `mgi`: Multigrid I-cycles

#### Generated Files
- `master-mac` - Compiled executable
- `slave.o` - Filtered object file with coupling symbols
- `vort-master.mp4` - Master domain vorticity visualization
- `vort-slave.mp4` - Slave domain vorticity visualization
- `*.ppm` - Individual animation frames (if ffmpeg fails)

## 🔍 Troubleshooting

### Docker Issues

#### Container Not Found
```bash
# Check container status
docker ps -a

# If stopped, restart
./docker_shell.sh start

# If doesn't exist, rebuild
./run_docker.sh
```

#### Docker Daemon Not Running
**Symptoms:** `Cannot connect to the Docker daemon`
**Solution:**
1. Open Docker Desktop application
2. Wait for Docker to start (whale icon in menu bar)
3. Try command again

#### Permission Issues
```bash
# If container permissions are wrong
docker exec basilisk-work chown -R $(id -u):$(id -g) /host
```

### Compilation Issues

#### `qcc: command not found`
**Cause:** Basilisk environment not set up
**Solution:**
```bash
# Inside container
export BASILISK=/workspace/basilisk-C/basilisk-source/src
export PATH=$PATH:$BASILISK

# Or restart container (auto-sets environment)
exit
./docker_shell.sh
```

#### `objcopy: command not found`
**Cause:** Not in Linux container (running on Mac)
**Solution:**
- Ensure you're inside Docker container (`./docker_shell.sh`)
- Use `/host` directory, not Mac terminal
- GNU binutils should be installed in container

#### Symbol Conflicts
**Error:** `multiple definition of symbols`
**Cause:** Running compilation on macOS without symbol filtering
**Solution:**
- Must use Docker environment
- GNU `objcopy` is required for symbol filtering
- Never compile master-slave coupling directly on macOS

#### Build Fails in Container
```bash
# Clean build
cd /host
rm -f slave.o master-mac *.o
./compile.sh
```

### Simulation Issues

#### No Output Files Generated
**Possible Causes:**
- Simulation crashed early
- Insufficient runtime for file generation
- Directory permission issues

**Solutions:**
```bash
# Check for error messages
./master-mac 2>&1 | tee simulation.log

# Run with timeout to see partial results
timeout 30s ./master-mac

# Check file permissions
ls -la /host/
```

#### Movies Not Generated
**Symptoms:** Simulation runs but no `.mp4` files created
**Causes:**
- `ffmpeg` not available (should be in container)
- `ppm2mp4` function issues
- Simulation stopped before movie generation

**Solutions:**
```bash
# Check if ppm files exist
ls -la *.ppm

# Manual movie creation (if ppm files exist)
ffmpeg -r 25 -i vort-slave-%04d.ppm vort-slave.mp4
```

#### Simulation Hangs
```bash
# Run with timeout
timeout 120s ./master-mac

# Check system resources
docker stats basilisk-work
```

### File Management Issues

#### Results Not Visible on Mac
**Cause:** Files created in container copy, not mounted directory
**Solution:**
- Always use `/host` directory in container
- Use `./copy_from_docker.sh` for organized results

#### Mount Issues
```bash
# Verify mount is working
./docker_shell.sh
ls -la /host/           # Should show your Mac files
touch /host/test.txt
exit
ls -la test.txt         # Should exist on Mac
rm test.txt
```

## 🎨 Customization & Extension

### Modifying Simulation Parameters

#### Edit Source Files (Mac)
```bash
# Edit with any Mac editor
vim master-mac.c
code slave-mac.c

# Common parameters to modify:
# - Grid resolution: LEVEL
# - Domain size: L0
# - Simulation time: T
# - Reynolds number: Re
```

#### Recompile and Test
```bash
./docker_shell.sh
cd /host
./compile.sh
timeout 60s ./master-mac    # Test with timeout
```

### Adding Dependencies to Container

Edit `Dockerfile`:
```dockerfile
RUN apt-get update && apt-get install -y \
    your-package-here \
    another-package \
    && rm -rf /var/lib/apt/lists/*
```

Rebuild:
```bash
docker build -t basilisk-ubuntu .
```

### Extending Scripts

All scripts follow these patterns:
- **Error handling:** Check preconditions and provide helpful messages
- **Help text:** Include usage examples and descriptions
- **Modularity:** Functions for reusable functionality
- **Logging:** Clear output for debugging

Example extension:
```bash
# Add to docker_shell.sh
debug)
    echo "Starting debug session..."
    docker exec -it basilisk-work gdb /host/master-mac
    ;;
```

### Performance Tuning

#### Container Resources
```bash
# Check current usage
docker stats basilisk-work

# Modify Docker Desktop settings:
# - Increase memory limit (default 2GB may be insufficient)
# - Increase CPU cores for parallel compilation
```

#### Simulation Performance
```bash
# Profile simulation
time ./master-mac

# Memory usage
docker exec basilisk-work ps aux | grep master-mac
```

## 📚 References & Resources

### Basilisk Documentation
- [Basilisk Official Website](http://basilisk.fr)
- [Basilisk Tutorial](http://basilisk.fr/Tutorial)
- [Master-Slave Coupling Example](http://basilisk.fr/src/examples/master.c)
- [Basilisk Google Group](https://groups.google.com/d/forum/basilisk-fr)

### Docker Resources
- [Docker Documentation](https://docs.docker.com/)
- [Docker Desktop for Mac](https://docs.docker.com/desktop/mac/)

### Development Tools
- [GNU Binutils Documentation](https://sourceware.org/binutils/docs/)
- [objcopy Manual](https://linux.die.net/man/1/objcopy)

## 💡 Pro Tips

### Development Efficiency
- **Keep container running:** Use `./docker_shell.sh` multiple times instead of restarting
- **Fast iteration:** Use `cd /host && ./compile.sh && timeout 30s ./master-mac` for quick tests
- **File organization:** Use `./copy_from_docker.sh` for organized output management
- **Editor integration:** Edit on Mac with full IDE support, compile in container

### Performance Optimization
- **Docker resources:** Increase memory/CPU limits in Docker Desktop
- **Parallel builds:** Container uses multiple cores automatically
- **File caching:** Docker layer caching speeds up rebuilds
- **Result management:** Copy only needed files to save disk space

### Debugging Strategies
- **Incremental testing:** Test compilation before full simulation runs
- **Timeout usage:** Use `timeout` command to limit run time during testing
- **Log analysis:** Save simulation output with `./master-mac 2>&1 | tee log.txt`
- **Container inspection:** Use `docker exec -it basilisk-work bash` for debugging

### Best Practices
- **Version control:** Keep source files in git, ignore build artifacts
- **Backup results:** Important simulation data should be copied and archived
- **Clean builds:** Regularly clean artifacts with `rm -f *.o master-mac`
- **Documentation:** Comment your parameter changes in source files

---

*This unified guide replaces CLEAN-SETUP.md, README-docker.md, and WORKFLOW.md*

*Last updated: August 2025*