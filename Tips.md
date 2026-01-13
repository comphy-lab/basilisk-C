# Tips and Troubleshooting

This document contains tips and solutions for common issues encountered when working with the Basilisk-C repository.

## Installing Darcs

Darcs is not available in the Ubuntu package repositories for newer Ubuntu versions. Instead, you can install it using Cabal:

```shell
# First install Cabal
sudo apt-get update
sudo apt-get install -y ghc cabal-install

# Then install Darcs
cabal update
cabal install darcs
```

If you encounter dependency conflicts, try installing a specific version of Darcs:

```shell
# Try different versions until one works
for version in 2.16.5 2.14.2 2.12.5 2.10.2; do
  echo "Trying Darcs version $version..."
  if cabal install darcs-$version; then
    echo "Successfully installed Darcs $version"
    break
  fi
done
```

After installation, make sure to add the Cabal binary directory to your PATH:

```shell
echo "$HOME/.cabal/bin" >> ~/.bashrc
source ~/.bashrc
```

## Installing OpenConnect VPN for Discoverer HPC on Ubuntu 24.04

If you encounter issues installing OpenConnect VPN on Ubuntu 24.04 (Noble), follow these steps:

1. Edit `/etc/apt/sources.list.d/ubuntu.sources` and replace "noble" with "jammy" at two instances:
   - Change `Suites: noble noble-updates noble-backports` to `Suites: jammy noble-updates noble-backports`
   - Change `Suites: noble-security` to `Suites: jammy-security`

2. Update and upgrade your system:
   ```shell
   sudo apt-get update
   sudo apt-get upgrade
   ```

3. Remove mailcap:
   ```shell
   sudo apt autoremove
   ```

4. Force install the required dependency:
   ```shell
   sudo apt install libtss2-mu0=3.2.0-1ubuntu1.1
   ```

5. Install OpenConnect:
   ```shell
   sudo apt-get install globalprotect-openconnect
   ```

6. **Important**: After installing OpenConnect, you may need to switch the sources back to "noble" to avoid dependency issues with other packages (like python3).

## GitHub Actions Workflow for Darcs Synchronization

We've set up a GitHub Actions workflow to automatically sync this repository with the Darcs repositories. The workflow:

1. Runs daily at midnight UTC
2. Can be triggered manually
3. Installs Darcs using Stack or Cabal as fallback
4. Pulls all changes from both Darcs repositories:
   - Syncs basilisk-source/ with http://basilisk.fr/basilisk
   - Syncs basilisk-wiki/ with http://basilisk.fr/wiki
5. Commits and pushes any changes to GitHub

The workflow file is located at `.github/workflows/sync-darcs-repositories.yml`.

## Installing Basilisk

We provide three installation scripts depending on your system and available tools:

### Option 1: Using Darcs (macOS with Homebrew)

```shell
./reset_install_basilisk.sh
```

This script:
1. Clones the latest Basilisk from http://basilisk.fr/basilisk using Darcs
2. Automatically downloads and applies all comphy-lab patches (including macOS compatibility)
3. Configures for your platform (`config.osx` or `config.gcc`)
4. Builds Basilisk

**Requirements**: `darcs`, `make`, `gcc`

### Option 2: Using Git (Recommended for most users)

```shell
./reset_install_basilisk-no-darcs.sh
```

This script:
1. Uses git sparse checkout to clone only `basilisk-source/` from our GitHub fork
2. Applies comphy-lab patches (macOS only)
3. Configures and builds Basilisk

**Requirements**: `git`, `make`, `gcc`, `gawk`

### Option 3: Using wget/tar (No version control needed)

```shell
./reset_install_basilisk-no-darcs-no-git.sh
```

This script:
1. Downloads Basilisk directly from http://basilisk.fr/basilisk/basilisk.tar.gz
2. Applies comphy-lab patches (macOS only)
3. Configures and builds Basilisk

**Requirements**: `wget`, `tar`, `make`, `gcc`, `gawk`

### Clean Reinstall

For any script, use the `--hard` flag to remove existing installation and start fresh:

```shell
./reset_install_basilisk.sh --hard
# or
./reset_install_basilisk-no-darcs.sh --hard
# or
./reset_install_basilisk-no-darcs-no-git.sh --hard
```

### Manual Installation

If you prefer to install manually:

```shell
# Clone Basilisk
darcs clone http://basilisk.fr/basilisk
cd basilisk

# Apply the macOS compatibility patch (macOS only)
curl -O https://raw.githubusercontent.com/comphy-lab/basilisk-C/main/patches/2025-11-03-macos-mman-compatibility.patch
patch -p1 < 2025-11-03-macos-mman-compatibility.patch

# Configure and build
cd src
ln -s config.osx config  # or config.gcc for Linux
make -k
make
```

### About the Patch

The `2025-11-03-macos-mman-compatibility.patch` fixes compilation issues on macOS by:
- Defining missing memory mapping constants (`MAP_ANON`, `MAP_ANONYMOUS`)
- Defining missing memory advice constants (`POSIX_MADV_DONTNEED`, `MADV_DONTNEED`)
- Explicitly declaring the `madvise()` function

**Credit**: Thanks to Peter Croxford for identifying this issue.

## Repository Structure

Our fork is organized with two main Darcs repositories:

- **basilisk-source/**: Contains the source code for Basilisk C
- **basilisk-wiki/**: Contains the documentation and wiki content

Each directory is synchronized daily with the corresponding Darcs repository from basilisk.fr.

## Reporting Issues

When reporting issues with Basilisk-C, please use our issue templates:

1. **Bug Reports**: For software bugs and unexpected behavior
2. **Installation Issues**: For problems encountered during installation
3. **Feature Requests**: For suggesting new features or improvements to build on top of Basilisk-C
4. **Generic Issues**: For other types of issues that don't fit the above categories

These templates help us gather all the necessary information to diagnose and fix issues quickly. When creating a new issue, you'll be prompted to choose the appropriate template.

## Using bview with a Local Client

By default, `bview` uses the remote JavaScript client hosted at basilisk.fr for 3D visualization. For offline development or faster loading, you can use a local client instead.

### Setup

1. Clone and serve the local client:
   ```shell
   git clone https://github.com/comphy-lab/bview-local-client.git
   cd bview-local-client
   # Start a local HTTP server (e.g., using Python)
   python3 -m http.server 8000
   ```

2. Run bview with the `--local` flag:
   ```shell
   bview2D --local dump      # Uses localhost:8000
   bview3D --local=3000 dump # Uses localhost:3000
   ```

### How It Works

The visualization URL has two parts:
- **HTTP server**: Serves the JavaScript 3D client (static files)
- **WebSocket**: Connects back to your local bview process

The `--local` flag only changes where the JS client is served from—the WebSocket always connects to your local bview process regardless.

For more details, see the [bview-local-client repository](https://github.com/comphy-lab/bview-local-client).
