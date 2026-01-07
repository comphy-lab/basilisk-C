#!/bin/bash
# Unified Basilisk installation script with multiple installation modes
# Tested on macOS and Linux. Report issues at https://github.com/comphy-lab/basilisk-C/issues
# Based on https://basilisk.fr/src/INSTALL

set -e

# ============================================================================
# Configuration
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$SCRIPT_DIR"
PROJECT_CONFIG="$REPO_ROOT/.project_config"
BASILISK_DIR="$REPO_ROOT/basilisk"
BASILISK_SRC_DIR="$BASILISK_DIR/src"
PATCHES_DIR="$REPO_ROOT/patches"

# GitHub URLs for patches
PATCHES_API_URL="https://api.github.com/repos/comphy-lab/basilisk-C/contents/patches"
PATCHES_RAW_URL="https://raw.githubusercontent.com/comphy-lab/basilisk-C/main/patches"

# ============================================================================
# Parse command line arguments
# ============================================================================

HARD_RESET=false
LOCAL_BVIEW=false
SHOW_HELP=false
MODE=""

for arg in "$@"; do
    case "$arg" in
        --hard)
            HARD_RESET=true
            ;;
        --local-bview)
            LOCAL_BVIEW=true
            ;;
        --help|-h)
            SHOW_HELP=true
            ;;
        --mode=*)
            MODE="${arg#*=}"
            ;;
    esac
done

# ============================================================================
# Help function
# ============================================================================

show_help() {
    cat << 'EOF'
Basilisk Installation Script
=============================

Usage: ./reset_install_requirements.sh [OPTIONS]

Options:
  --help, -h      Show this help message
  --hard          Force reinstall (removes existing basilisk directory)
  --local-bview   Apply the local-bview patch for localhost JavaScript client
  --mode=N        Select installation mode non-interactively (1-3)

Installation Modes:
  1) default       - Use darcs to clone basilisk, fetch patches from GitHub
                     Best for: Standard development workflow
                     Requires: darcs, make, gcc

  2) remote-fr     - Download tarball from basilisk.fr, fetch patches from GitHub
                     Best for: Systems without darcs installed
                     Requires: wget, tar, curl, gawk, make, gcc

  3) remote-comphy - Clone from comphy-lab GitHub fork (pre-patched)
                     Best for: Quick setup without darcs
                     Requires: git, curl, gawk, make, gcc

Examples:
  ./reset_install_requirements.sh                    # Interactive mode selection
  ./reset_install_requirements.sh --mode=1           # Use default mode
  ./reset_install_requirements.sh --mode=2 --hard    # Reinstall using wget
  ./reset_install_requirements.sh --mode=1 --local-bview  # Include local-bview patch

Note: For testing local patches before pushing to GitHub, use test-patch-local.sh instead.

For more information, visit: https://github.com/comphy-lab/basilisk-C
EOF
}

# ============================================================================
# Utility functions
# ============================================================================

print_green() {
    printf "\033[0;32m%s\033[0m\n" "$1"
}

print_red() {
    printf "\033[0;31m%s\033[0m\n" "$1"
}

print_yellow() {
    printf "\033[0;33m%s\033[0m\n" "$1"
}

print_cyan() {
    printf "\033[0;36m%s\033[0m\n" "$1"
}

# ============================================================================
# Prerequisite checking functions
# ============================================================================

check_tool() {
    local tool="$1"
    if command -v "$tool" > /dev/null 2>&1; then
        print_green "✓ $tool is installed"
        return 0
    else
        return 1
    fi
}

show_install_instructions() {
    local tool="$1"

    if [[ "$OSTYPE" == "darwin"* ]]; then
        case "$tool" in
            make|gcc)
                echo "  xcode-select --install"
                ;;
            darcs)
                echo "  brew install darcs"
                ;;
            wget)
                echo "  brew install wget"
                ;;
            gawk)
                echo "  brew install gawk"
                ;;
            git)
                echo "  xcode-select --install"
                ;;
            curl)
                echo "  brew install curl"
                ;;
            tar)
                echo "  (tar should be pre-installed on macOS)"
                ;;
        esac
    else
        case "$tool" in
            make|gcc)
                echo "  sudo apt install build-essential"
                ;;
            darcs)
                echo "  Visit https://darcs.net/ for installation instructions"
                ;;
            wget|gawk|git|curl|tar)
                echo "  sudo apt install $tool"
                ;;
        esac
    fi
}

check_prerequisites() {
    local mode="$1"
    local missing_tools=()

    echo "Checking prerequisites for mode $mode..."
    echo ""

    # Common prerequisites
    check_tool "make" || missing_tools+=("make")
    check_tool "gcc" || missing_tools+=("gcc")

    # Mode-specific prerequisites
    case "$mode" in
        1)    # darcs mode
            check_tool "darcs" || missing_tools+=("darcs")
            ;;
        2)    # wget mode
            check_tool "wget" || missing_tools+=("wget")
            check_tool "tar" || missing_tools+=("tar")
            check_tool "curl" || missing_tools+=("curl")
            check_tool "gawk" || missing_tools+=("gawk")
            ;;
        3)    # git mode
            check_tool "git" || missing_tools+=("git")
            check_tool "curl" || missing_tools+=("curl")
            check_tool "gawk" || missing_tools+=("gawk")
            ;;
    esac

    echo ""

    if [[ ${#missing_tools[@]} -gt 0 ]]; then
        print_red "Error: Missing required tools: ${missing_tools[*]}"
        echo ""
        echo "Installation instructions:"
        for tool in "${missing_tools[@]}"; do
            echo "  $tool:"
            show_install_instructions "$tool"
        done
        echo ""
        echo "Please install the missing tools and try again."
        exit 1
    else
        print_green "✅ All prerequisites are satisfied!"
        echo ""
    fi
}

# ============================================================================
# Mode selection
# ============================================================================

select_mode() {
    echo ""
    print_cyan "Basilisk Installation Options:"
    echo "  1) default       - darcs clone + GitHub patches (recommended)"
    echo "  2) remote-fr     - wget tarball + GitHub patches (no darcs needed)"
    echo "  3) remote-comphy - git clone from comphy-lab fork (no darcs needed)"
    echo ""
    printf "Select mode [1-3, default=1]: "
    read -r selection

    # Default to mode 1 if empty
    if [[ -z "$selection" ]]; then
        selection="1"
    fi

    # Validate selection
    case "$selection" in
        1|2|3)
            MODE="$selection"
            ;;
        *)
            print_red "Invalid selection: $selection"
            exit 1
            ;;
    esac

    echo ""
    print_cyan "Selected mode: $MODE"
    echo ""
}

# ============================================================================
# Patch application functions
# ============================================================================

apply_patches_github() {
    local target_dir="$1"
    local apply_local_bview="${2:-false}"
    local patch_failed=false

    if [[ "$OSTYPE" != "darwin"* ]]; then
        # Patches are macOS-specific, skip on other platforms
        return 0
    fi

    print_cyan "Applying comphy-lab patches from GitHub..."

    # Create temp directory for patches
    mkdir -p "$target_dir/.patches_temp"

    # Get list of patch files (sorted by name for chronological order due to YYYY-MM-DD format)
    local PATCH_FILES
    PATCH_FILES=$(curl -s "$PATCHES_API_URL" | grep -o '"name": "[^"]*\.patch"' | sed 's/"name": "//;s/"//' | sort)

    if [[ -z "$PATCH_FILES" ]]; then
        print_yellow "Warning: No patches found or unable to fetch patch list"
    else
        # Download and apply each patch
        while read -r patch_file; do
            if [[ -n "$patch_file" ]]; then
                # Skip local-bview patch unless --local-bview flag was provided
                if [[ "$patch_file" == *"-local-bview.patch" ]] && [[ "$apply_local_bview" != "true" ]]; then
                    echo "  Skipping $patch_file (use --local-bview to apply)"
                    continue
                fi
                echo "  Downloading $patch_file..."
                if curl -s -f "$PATCHES_RAW_URL/$patch_file" -o "$target_dir/.patches_temp/$patch_file"; then
                    echo "  Applying $patch_file..."
                    if (cd "$target_dir" && patch -p1 < ".patches_temp/$patch_file"); then
                        print_green "  ✓ Successfully applied $patch_file"
                    else
                        print_red "  ✗ Failed to apply $patch_file"
                        patch_failed=true
                    fi
                else
                    print_red "  ✗ Failed to download $patch_file"
                    patch_failed=true
                fi
            fi
        done <<< "$PATCH_FILES"
    fi

    # Clean up
    rm -rf "$target_dir/.patches_temp"
    echo ""

    if [[ "$patch_failed" == true ]]; then
        return 1
    fi
    return 0
}

# ============================================================================
# Basilisk installation functions
# ============================================================================

install_basilisk_darcs() {
    print_cyan "Cloning basilisk using darcs..."

    if ! darcs clone https://basilisk.fr/basilisk "$BASILISK_DIR"; then
        print_red "Error: Failed to clone basilisk into $BASILISK_DIR"
        exit 1
    fi
}

install_basilisk_wget() {
    print_cyan "Downloading basilisk using wget..."

    if ! cd "$REPO_ROOT"; then
        print_red "Error: Failed to change directory to $REPO_ROOT"
        exit 1
    fi

    if ! wget https://basilisk.fr/basilisk/basilisk.tar.gz; then
        print_red "Error: Failed to download basilisk.tar.gz"
        exit 1
    fi

    print_cyan "Extracting basilisk.tar.gz..."
    if ! tar xzf basilisk.tar.gz; then
        print_red "Error: Failed to extract basilisk.tar.gz"
        exit 1
    fi

    # Clean up the tar file
    rm basilisk.tar.gz
}

install_basilisk_git() {
    local REPO_URL="https://github.com/comphy-lab/basilisk-C.git"
    local TEMP_DIR="$REPO_ROOT/basilisk-C-temp"

    print_cyan "Cloning basilisk-source from comphy-lab/basilisk-C (sparse checkout)..."

    # Use sparse checkout to only get basilisk-source directory
    if ! git clone --depth 1 --filter=blob:none --sparse "$REPO_URL" "$TEMP_DIR"; then
        print_red "Error: Failed to clone repository"
        echo "URL: $REPO_URL"
        exit 1
    fi

    if ! cd "$TEMP_DIR"; then
        print_red "Error: Failed to change directory to $TEMP_DIR"
        exit 1
    fi

    if ! git sparse-checkout set basilisk-source; then
        print_red "Error: Failed to set sparse checkout for basilisk-source"
        exit 1
    fi

    if ! cd "$REPO_ROOT"; then
        print_red "Error: Failed to change directory to $REPO_ROOT"
        exit 1
    fi

    # Move basilisk-source to basilisk
    print_cyan "Setting up basilisk directory..."
    if ! mv "$TEMP_DIR/basilisk-source" "$BASILISK_DIR"; then
        print_red "Error: Failed to move basilisk-source into $BASILISK_DIR"
        exit 1
    fi

    # Clean up temp clone
    rm -rf "$TEMP_DIR"
}

# ============================================================================
# Build function
# ============================================================================

build_basilisk() {
    if ! cd "$BASILISK_SRC_DIR"; then
        print_red "Error: Failed to change directory to $BASILISK_SRC_DIR"
        exit 1
    fi

    # Link appropriate config file
    if [[ "$OSTYPE" == "darwin"* ]]; then
        print_cyan "Using macOS configuration..."
        ln -s config.osx config
    else
        print_cyan "Using Linux configuration..."
        ln -s config.gcc config
    fi

    print_cyan "Building basilisk (first pass with -k to continue on errors)..."
    if ! make -k; then
        print_red "Error: make -k failed in $BASILISK_SRC_DIR"
        exit 1
    fi

    print_cyan "Building basilisk (final build)..."
    if ! make; then
        print_red "Error: make failed in $BASILISK_SRC_DIR"
        exit 1
    fi
}

# ============================================================================
# Environment setup
# ============================================================================

setup_environment() {
    if [[ ! -d "$BASILISK_SRC_DIR" ]]; then
        print_red "Error: Expected $BASILISK_SRC_DIR to exist before writing .project_config"
        exit 1
    fi

    printf "export BASILISK=%s\n" "$BASILISK_SRC_DIR" > "$PROJECT_CONFIG"
    # Prepend BASILISK to PATH so Basilisk tools (qcc, etc.) take precedence
    printf "export PATH=\$BASILISK:\$PATH\n" >> "$PROJECT_CONFIG"

    source "$PROJECT_CONFIG"
}

verify_installation() {
    echo ""
    print_cyan "Checking qcc installation..."

    if ! qcc --version > /dev/null 2>&1; then
        print_red "Error: qcc is not working properly."
        if [[ "$OSTYPE" == "darwin"* ]]; then
            echo "Please ensure you have Xcode Command Line Tools installed."
            echo "You can install them by running: xcode-select --install"
        else
            echo "Please ensure you have build-essential installed."
            echo "You can install it by running: sudo apt install build-essential"
        fi
        echo "For more details, visit: https://basilisk.fr/src/INSTALL"
        exit 1
    else
        print_green "✅ qcc is properly installed."
        qcc --version
    fi
}

# ============================================================================
# Main execution
# ============================================================================

# Show help if requested
if [[ "$SHOW_HELP" == true ]]; then
    show_help
    exit 0
fi

# Select mode interactively if not provided
if [[ -z "$MODE" ]]; then
    select_mode
fi

# Validate mode
case "$MODE" in
    1|2|3)
        ;;
    *)
        print_red "Invalid mode: $MODE (must be 1-3)"
        exit 1
        ;;
esac

# Check prerequisites for selected mode
check_prerequisites "$MODE"

# Remove project config always
rm -rf "$PROJECT_CONFIG"

# Install basilisk based on mode
if [[ "$HARD_RESET" == true ]] || [[ ! -d "$BASILISK_DIR" ]]; then
    print_cyan "Installing basilisk (mode $MODE)..."
    rm -rf "$BASILISK_DIR"

    case "$MODE" in
        1)  # default: darcs + GitHub patches
            install_basilisk_darcs
            if ! apply_patches_github "$BASILISK_DIR" "$LOCAL_BVIEW"; then
                print_red "Error: Failed to apply patches"
                exit 1
            fi
            ;;
        2)  # remote-fr: wget + GitHub patches
            install_basilisk_wget
            if ! apply_patches_github "$BASILISK_DIR" "$LOCAL_BVIEW"; then
                print_red "Error: Failed to apply patches"
                exit 1
            fi
            ;;
        3)  # remote-comphy: git sparse checkout + GitHub patches
            install_basilisk_git
            if ! apply_patches_github "$BASILISK_DIR" "$LOCAL_BVIEW"; then
                print_red "Error: Failed to apply patches"
                exit 1
            fi
            ;;
    esac

    build_basilisk
else
    print_cyan "Using existing basilisk installation..."
    if [[ ! -d "$BASILISK_SRC_DIR" ]]; then
        print_red "Error: Missing $BASILISK_SRC_DIR. Run with --hard to reinstall."
        exit 1
    fi
fi

# Setup environment and verify
setup_environment
verify_installation

echo ""
print_green "✅ Basilisk installation complete!"
echo "To use basilisk in your shell, run:"
echo "  source $PROJECT_CONFIG"
