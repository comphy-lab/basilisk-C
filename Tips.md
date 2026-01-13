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

We provide four installation scripts depending on your system and available tools:

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

### Option 4: Ref-Locked Install (Reproducible, HPC-friendly)

Use this when you want to lock Basilisk + patches to a GitHub **Release tag** in `comphy-lab/basilisk-C`.

```shell
./reset_install_basilisk-ref-locked.sh --ref=v2026-01-13 --hard
```

Equivalent (using the unified installer):

```shell
./reset_install_basilisk.sh --mode=4 --ref=v2026-01-13 --hard
```

This script:
1. Downloads the OS-specific tarball attached to the Release tag (`basilisk-mac.tar.gz` or `basilisk-linux.tar.gz`)
2. Builds Basilisk and writes a lock stamp to `basilisk/.comphy-lock`

**Requirements**: `curl`, `tar`, `make`, `gcc`, `gawk` (plus `patch` only if using `--local-bview`)

### Creating a Ref-Locked Release (Maintainers)

This repo includes `release-comphy-tag.sh` to create a reproducible, ref-locked Release tag and publish OS-specific, pre-patched Basilisk source tarballs.

```shell
./release-comphy-tag.sh
# or
./release-comphy-tag.sh --tag=v2026-01-13
```

This script:
1. Pulls the latest upstream Basilisk snapshot into `basilisk-source/` using `darcs`
2. Runs install tests (`reset_install_basilisk.sh --mode=1 --hard`):
   - macOS: tests macOS natively + Linux in Docker (`darcs-test` image)
   - Linux: tests Linux natively
3. Builds and uploads GitHub Release assets:
   - `basilisk-mac.tar.gz` (+ `.sha256`): upstream snapshot + macOS + common patches (**local-bview not applied**)
   - `basilisk-linux.tar.gz` (+ `.sha256`): upstream snapshot + common patches (**local-bview not applied**, macOS patches excluded)

### Clean Reinstall

For any script, use the `--hard` flag to remove existing installation and start fresh:

```shell
./reset_install_basilisk.sh --hard
# or
./reset_install_basilisk-no-darcs.sh --hard
# or
./reset_install_basilisk-no-darcs-no-git.sh --hard
# or (ref-locked)
./reset_install_basilisk-ref-locked.sh --ref=v2026-01-13 --hard
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

## Patches

Our fork maintains a set of patches in the `patches/` directory that fix issues or add features not yet available in upstream Basilisk. These patches are automatically applied by the installation scripts. For detailed information about each patch, see [`patches/README.md`](patches/README.md).

### Available Patches

| Patch | Platform | Description |
|-------|----------|-------------|
| `2025-11-03-macos-mman-compatibility.patch` | macOS | Fixes memory mapping compilation errors |
| `2026-01-06-local-bview.patch` | All | Adds `--local` flag to bview for offline visualization |
| `2026-01-13-mpi-tree-dump-header-fix.patch` | All | Fixes uninitialized header in MPI tree dump |

### Applying Patches Manually

If you need to apply patches manually to an existing Basilisk installation:

```shell
cd /path/to/basilisk
curl -O https://raw.githubusercontent.com/comphy-lab/basilisk-C/main/patches/PATCH_NAME.patch
patch -p1 < PATCH_NAME.patch
cd src && make clean && make
```

## Repository Structure

Our fork is organized as follows:

- **basilisk-source/**: Contains the source code for Basilisk C (synced from Darcs)
- **basilisk-wiki/**: Contains the documentation and wiki content (synced from Darcs)
- **patches/**: comphy-lab patches for bug fixes and enhancements
- **bugs/**: Test cases and minimal reproductions for tracking issues

The `basilisk-source/` and `basilisk-wiki/` directories are synchronized daily with the corresponding Darcs repositories from basilisk.fr.

## Reporting Issues

When reporting issues with Basilisk-C, please use our issue templates:

1. **Bug Reports**: For software bugs and unexpected behavior
2. **Installation Issues**: For problems encountered during installation
3. **Feature Requests**: For suggesting new features or improvements to build on top of Basilisk-C
4. **Generic Issues**: For other types of issues that don't fit the above categories

These templates help us gather all the necessary information to diagnose and fix issues quickly. When creating a new issue, you'll be prompted to choose the appropriate template.

## Using bview with a Local Client

By default, `bview2D`/`bview3D` prints a URL which uses the remote JavaScript client (served from `basilisk.ida.upmc.fr`). For offline use or faster loading, you can serve the client locally using [bview-local-client](https://github.com/comphy-lab/bview-local-client).

### Setup

1. Clone and serve the local client:
   ```shell
   git clone https://github.com/comphy-lab/bview-local-client.git
   cd bview-local-client
   ./deploy.sh  # serves http://localhost:8000/three.js/editor/index.html
   ```

### Option A: No Basilisk patch (manual URL)

This works with upstream Basilisk as-is: run `bview` normally, then use the local client URL and keep the `?ws://...` part from the printed output.

```shell
bview2D dump
# Example output:
#   http://basilisk.ida.upmc.fr/three.js/editor/index.html?ws://your-hostname:7101

# Open the local client and append the websocket query:
#   http://localhost:8000/three.js/editor/index.html?ws://your-hostname:7101
```

### Option B: Optional convenience patch (prints localhost URL)

`patches/2026-01-06-local-bview.patch` is an optional, aesthetic/convenience patch: it only changes the JavaScript client base URL that `bview` prints (and what `display.html` redirects to). It does **not** change the websocket connection (which always connects back to your local `bview` process).

To enable it, install Basilisk with `--local-bview` (not applied by default), then run `bview` with the `--local` flag:
   ```shell
   # installer examples
   ./reset_install_basilisk.sh --mode=1 --local-bview --hard
   ./reset_install_basilisk.sh --mode=4 --ref=v2026-01-13 --local-bview --hard

   # usage
   bview2D --local dump      # Uses localhost:8000
   bview3D --local=3000 dump # Uses localhost:3000
   ```
   **Note:** `--local` hardcodes `localhost`. If you want to open the visualization from another machine, use Option A and replace the base URL with the correct IP/hostname.

### How It Works

The visualization URL has two parts:
- **HTTP server**: Serves the JavaScript 3D client (static files)
- **WebSocket**: Connects back to your local bview process

The local-bview patch only changes where the JS client is served from—the websocket always connects to your local bview process regardless.

**Note:** Our GitHub Release tarballs (`basilisk-mac.tar.gz` / `basilisk-linux.tar.gz`) intentionally exclude `2026-01-06-local-bview.patch`. If you want `bview --local` with a ref-locked install, pass `--local-bview` so the installer downloads and applies the patch for that same `--ref`.

For more details, see the [bview-local-client repository](https://github.com/comphy-lab/bview-local-client).
