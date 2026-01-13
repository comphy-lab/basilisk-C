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
REF=""

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
        --ref=*)
            REF="${arg#*=}"
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

Usage: ./reset_install_basilisk.sh [OPTIONS]

Options:
  --help, -h      Show this help message
  --hard          Force reinstall (removes existing basilisk directory)
  --local-bview   Apply optional local-bview convenience patch (enables `bview --local` URL output for bview-local-client)
  --mode=N        Select installation mode non-interactively (1-4)
  --ref=REF       Pinned tag (release) for mode 4 only

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

  4) ref-locked    - Install from OS-specific tarball attached to a GitHub Release tag
                     Best for: Reproducible installs, HPC (no darcs/git)
                     Requires: curl, tar, gawk, make, gcc

Examples:
  ./reset_install_basilisk.sh                    # Interactive mode selection
  ./reset_install_basilisk.sh --mode=1           # Use default mode
  ./reset_install_basilisk.sh --mode=2 --hard    # Reinstall using wget
  ./reset_install_basilisk.sh --mode=1 --local-bview  # Include local-bview patch
  ./reset_install_basilisk.sh --mode=4 --ref=v2026-01-13 --hard  # Ref-locked install (from release assets)
  ./reset_install_basilisk.sh --mode=4 --ref=v2026-01-13 --local-bview --hard  # Ref-locked + local-bview patch

Notes:
  - The local-bview patch is optional and not applied by default.
  - GitHub Release tarballs intentionally exclude it; `--local-bview` downloads and applies the patch for the same `--ref`.

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
            patch)
                echo "  (patch should be pre-installed on macOS)"
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
            patch)
                echo "  sudo apt install patch"
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
            check_tool "patch" || missing_tools+=("patch")
            ;;
        2)    # wget mode
            check_tool "wget" || missing_tools+=("wget")
            check_tool "tar" || missing_tools+=("tar")
            check_tool "curl" || missing_tools+=("curl")
            check_tool "gawk" || missing_tools+=("gawk")
            check_tool "patch" || missing_tools+=("patch")
            ;;
        3)    # git mode
            check_tool "git" || missing_tools+=("git")
            check_tool "curl" || missing_tools+=("curl")
            check_tool "gawk" || missing_tools+=("gawk")
            check_tool "patch" || missing_tools+=("patch")
            ;;
        4)    # ref-locked mode
            check_tool "curl" || missing_tools+=("curl")
            check_tool "tar" || missing_tools+=("tar")
            check_tool "gawk" || missing_tools+=("gawk")
            if [[ "$LOCAL_BVIEW" == true ]]; then
                check_tool "patch" || missing_tools+=("patch")
            fi
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
    echo "  4) ref-locked    - pinned GitHub Release tag (reproducible)"
    echo ""
    printf "Select mode [1-4, default=1]: "
    read -r selection

    # Default to mode 1 if empty
    if [[ -z "$selection" ]]; then
        selection="1"
    fi

    # Validate selection
    case "$selection" in
        1|2|3|4)
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

    if [[ "$MODE" == "4" ]] && [[ -z "${REF}" ]]; then
        printf "Enter GitHub Release tag for mode 4 (e.g. v2026-01-13): "
        read -r REF
        if [[ -z "${REF}" ]]; then
            print_red "Error: --ref is required for mode 4"
            exit 1
        fi
    fi
}

read_lock_ref() {
    local lock_file="$1"
    local lock_ref=""

    if [[ ! -f "$lock_file" ]]; then
        echo ""
        return 0
    fi

    while IFS= read -r line; do
        case "$line" in
            ref=*)
                lock_ref="${line#ref=}"
                ;;
        esac
    done < "$lock_file"

    echo "$lock_ref"
}

read_lock_os() {
    local lock_file="$1"
    local lock_os=""

    if [[ ! -f "$lock_file" ]]; then
        echo ""
        return 0
    fi

    while IFS= read -r line; do
        case "$line" in
            os=*)
                lock_os="${line#os=}"
                ;;
        esac
    done < "$lock_file"

    echo "$lock_os"
}

write_lock_stamp() {
    local lock_file="$1"
    local ref="$2"
    local patches_dir="$3"
    local apply_local_bview="$4"

    local patch_files=()
    local applied_patches=()
    local skipped_patches=()

    patch_files=($(ls "$patches_dir"/*.patch 2>/dev/null | sort))
    for patch_file in "${patch_files[@]}"; do
        local patch_name
        patch_name=$(basename "$patch_file")

        if [[ "$patch_name" == *"-local-bview.patch" ]] && [[ "$apply_local_bview" != "true" ]]; then
            skipped_patches+=("$patch_name")
            continue
        fi
        if [[ "$patch_name" == *"-macos-"* ]] && [[ "$OSTYPE" != "darwin"* ]]; then
            skipped_patches+=("$patch_name")
            continue
        fi
        applied_patches+=("$patch_name")
    done

    {
        printf "ref=%s\n" "$ref"
        printf "os=%s\n" "$([[ "$OSTYPE" == "darwin"* ]] && echo "darwin" || echo "linux")"
        printf "created_utc=%s\n" "$(date -u +"%Y-%m-%dT%H:%M:%SZ" 2>/dev/null || date)"
        printf "local_bview=%s\n" "$apply_local_bview"
        printf "patches_applied=%s\n" "${applied_patches[*]}"
        printf "patches_skipped=%s\n" "${skipped_patches[*]}"
    } > "$lock_file"
}

apply_patches_from_dir() {
    local target_dir="$1"
    local patches_dir="$2"
    local apply_local_bview="${3:-false}"
    local patch_failed=false

    print_cyan "Applying comphy-lab patches (from pinned ref)..."

    if [[ ! -d "$patches_dir" ]]; then
        print_red "Error: Missing patches directory: $patches_dir"
        return 1
    fi

    local patch_files
    patch_files=($(ls "$patches_dir"/*.patch 2>/dev/null | sort))

    if [[ ${#patch_files[@]} -eq 0 ]]; then
        print_yellow "Warning: No patches found in $patches_dir"
        return 0
    fi

    for patch_file in "${patch_files[@]}"; do
        local patch_name
        patch_name=$(basename "$patch_file")

        if [[ "$patch_name" == *"-local-bview.patch" ]] && [[ "$apply_local_bview" != "true" ]]; then
            echo "  Skipping $patch_name (use --local-bview to apply)"
            continue
        fi
        if [[ "$patch_name" == *"-macos-"* ]] && [[ "$OSTYPE" != "darwin"* ]]; then
            echo "  Skipping $patch_name (macOS-specific patch)"
            continue
        fi

        echo "  Applying $patch_name..."
        if (cd "$target_dir" && patch -p1 < "$patch_file"); then
            print_green "  ✓ Successfully applied $patch_name"
        else
            print_red "  ✗ Failed to apply $patch_name"
            patch_failed=true
        fi
    done

    echo ""

    if [[ "$patch_failed" == true ]]; then
        return 1
    fi
    return 0
}

# ============================================================================
# Patch application functions
# ============================================================================

# Checks for local patches/ directory first, falls back to GitHub if not available
apply_patches_github() {
    local target_dir="$1"
    local apply_local_bview="${2:-false}"
    local patch_failed=false

    print_cyan "Applying comphy-lab patches..."

    # Check for local patches directory first (useful for HPC/offline environments)
    if [[ -d "$PATCHES_DIR" ]]; then
        print_cyan "Using local patches from: $PATCHES_DIR"

        # Get list of patch files sorted by name (YYYY-MM-DD format ensures chronological order)
        local patch_files
        patch_files=($(ls "$PATCHES_DIR"/*.patch 2>/dev/null | sort))

        if [[ ${#patch_files[@]} -eq 0 ]]; then
            print_yellow "Warning: No patches found in $PATCHES_DIR"
        else
            for patch_file in "${patch_files[@]}"; do
                local patch_name=$(basename "$patch_file")

                # Skip local-bview patch unless --local-bview flag was provided
                if [[ "$patch_name" == *"-local-bview.patch" ]] && [[ "$apply_local_bview" != "true" ]]; then
                    echo "  Skipping $patch_name (use --local-bview to apply)"
                    continue
                fi
                # Skip macOS-specific patches on non-macOS systems
                if [[ "$patch_name" == *"-macos-"* ]] && [[ "$OSTYPE" != "darwin"* ]]; then
                    echo "  Skipping $patch_name (macOS-specific patch)"
                    continue
                fi

                echo "  Applying $patch_name..."
                if (cd "$target_dir" && patch -p1 < "$patch_file"); then
                    print_green "  ✓ Successfully applied $patch_name"
                else
                    print_red "  ✗ Failed to apply $patch_name"
                    patch_failed=true
                fi
            done
        fi
    else
        # Fall back to downloading patches from GitHub (requires internet)
        print_cyan "No local patches/ directory found, downloading from GitHub..."
        print_yellow "Note: Internet connection required. For offline use, include patches/ directory."

        # Create temp directory for patches
        mkdir -p "$target_dir/.patches_temp"

        # Get list of patch files (sorted by name for chronological order due to YYYY-MM-DD format)
        local PATCH_FILES
        PATCH_FILES=$(curl -s "$PATCHES_API_URL" | grep -o '"name": "[^"]*\.patch"' | sed 's/"name": "//;s/"//' | sort)

        if [[ -z "$PATCH_FILES" ]]; then
            print_yellow "Warning: No patches found or unable to fetch patch list (check internet connection)"
        else
            # Download and apply each patch
            while read -r patch_file; do
                if [[ -n "$patch_file" ]]; then
                    # Skip local-bview patch unless --local-bview flag was provided
                    if [[ "$patch_file" == *"-local-bview.patch" ]] && [[ "$apply_local_bview" != "true" ]]; then
                        echo "  Skipping $patch_file (use --local-bview to apply)"
                        continue
                    fi
                    # Skip macOS-specific patches on non-macOS systems
                    if [[ "$patch_file" == *"-macos-"* ]] && [[ "$OSTYPE" != "darwin"* ]]; then
                        echo "  Skipping $patch_file (macOS-specific patch)"
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
    fi

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

install_basilisk_ref_locked() {
    local ref="$1"
    local lock_file="$2"

    local asset_name
    local download_url
    local temp_dir

    if [[ "$OSTYPE" == "darwin"* ]]; then
        asset_name="basilisk-mac.tar.gz"
    else
        asset_name="basilisk-linux.tar.gz"
    fi

    download_url="https://github.com/comphy-lab/basilisk-C/releases/download/${ref}/${asset_name}"

    print_cyan "Installing Basilisk from ref-locked release: $ref"
    print_cyan "Downloading: $download_url"

    temp_dir="$(mktemp -d 2>/dev/null || mktemp -d -t basilisk-C)"
    if [[ -z "$temp_dir" ]] || [[ ! -d "$temp_dir" ]]; then
        print_red "Error: Failed to create temp directory"
        exit 1
    fi

    if ! curl -s -f -L "$download_url" -o "$temp_dir/$asset_name"; then
        print_red "Error: Failed to download release asset for ref '$ref'"
        print_red "Expected GitHub Release asset: $asset_name"
        print_red "Hint: Create the release using ./release-comphy-tag.sh"
        rm -rf "$temp_dir"
        exit 1
    fi

    if ! tar xzf "$temp_dir/$asset_name" -C "$temp_dir"; then
        print_red "Error: Failed to extract downloaded tarball"
        rm -rf "$temp_dir"
        exit 1
    fi

    if [[ ! -d "$temp_dir/basilisk/src" ]]; then
        print_red "Error: Tarball is missing expected basilisk/src directory"
        print_red "Expected tarball layout: basilisk/ and patches/ at top-level"
        rm -rf "$temp_dir"
        exit 1
    fi

    print_cyan "Setting up basilisk directory..."
    if ! mv "$temp_dir/basilisk" "$BASILISK_DIR"; then
        print_red "Error: Failed to move basilisk into $BASILISK_DIR"
        rm -rf "$temp_dir"
        exit 1
    fi

    # Always use .comphy-patches for ref-locked installs (don't pollute repo root)
    local patches_dest="$BASILISK_DIR/.comphy-patches"
    if [[ -d "$temp_dir/patches" ]]; then
        if ! mv "$temp_dir/patches" "$patches_dest"; then
            print_red "Error: Failed to move patches into $patches_dest"
            rm -rf "$temp_dir"
            exit 1
        fi
    else
        print_red "Error: Tarball is missing expected patches/ directory"
        rm -rf "$temp_dir"
        exit 1
    fi

    if [[ "$LOCAL_BVIEW" == true ]]; then
        print_cyan "Applying local-bview patch (from pinned ref)..."

        local patches_api_url_ref="${PATCHES_API_URL}?ref=${ref}"
        local raw_base_url="https://raw.githubusercontent.com/comphy-lab/basilisk-C/${ref}/patches"
        local patch_file
        patch_file="$(curl -s "$patches_api_url_ref" | grep -o '\"name\": \"[^\"]*\.patch\"' | sed 's/\"name\": \"//;s/\"//' | grep -E -- '-local-bview\.patch$' | head -n 1)"

        if [[ -z "$patch_file" ]]; then
            print_red "Error: Could not find a *-local-bview.patch file at ref '$ref'"
            rm -rf "$temp_dir"
            exit 1
        fi

        if ! curl -s -f -L "$raw_base_url/$patch_file" -o "$patches_dest/$patch_file"; then
            print_red "Error: Failed to download $patch_file from pinned ref"
            rm -rf "$temp_dir"
            exit 1
        fi

        if ! (cd "$BASILISK_DIR" && patch -p1 < "$patches_dest/$patch_file"); then
            print_red "Error: Failed to apply $patch_file"
            rm -rf "$temp_dir"
            exit 1
        fi
    fi

    write_lock_stamp "$lock_file" "$ref" "$patches_dest" "$LOCAL_BVIEW"

    rm -rf "$temp_dir"
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

    {
        printf 'export BASILISK="%s"\n' "$BASILISK_SRC_DIR"
        cat <<'EOF'

if [ ! -d "$BASILISK" ]; then
  echo "Error: BASILISK directory not found: $BASILISK" >&2
  return 1 2>/dev/null || exit 1
fi

# Ensure the repo-installed Basilisk takes precedence, without duplicating PATH
_PATH=":$PATH:"
_PATH="${_PATH//:$BASILISK:/:}"
_PATH="${_PATH#:}"
_PATH="${_PATH%:}"
export PATH="$BASILISK${_PATH:+:$_PATH}"
unset _PATH

# Refresh command lookup cache (best-effort)
hash -r 2>/dev/null || true
rehash 2>/dev/null || true
EOF
    } > "$PROJECT_CONFIG"

    source "$PROJECT_CONFIG"
}

verify_installation() {
    echo ""
    print_cyan "Checking qcc installation..."

    local qcc_path
    qcc_path="$(command -v qcc 2>/dev/null || true)"

    if [[ -z "$qcc_path" ]]; then
        print_red "Error: qcc not found on PATH after sourcing .project_config"
        exit 1
    fi

    if [[ "$qcc_path" != /* ]]; then
        print_red "Error: qcc resolves to a non-path ('$qcc_path'). This is likely an alias/function."
        print_red "Remove/disable the alias/function and try again."
        exit 1
    fi

    if [[ "$qcc_path" != "$BASILISK/"* ]]; then
        print_red "Error: qcc resolves to '$qcc_path' but expected it under '$BASILISK/'."
        print_red "This usually means a global Basilisk install is taking precedence."
        print_red "Try: 'hash -r' (bash) or 'rehash' (zsh), then re-source .project_config."
        exit 1
    fi

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
    1|2|3|4)
        ;;
    *)
        print_red "Invalid mode: $MODE (must be 1-4)"
        exit 1
        ;;
esac

if [[ "$MODE" == "4" ]] && [[ -z "${REF}" ]]; then
    print_red "Error: --ref is required for mode 4"
    exit 1
fi

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
        4)  # ref-locked: pinned comphy-lab/basilisk-C ref
            install_basilisk_ref_locked "$REF" "$BASILISK_DIR/.comphy-lock"
            ;;
    esac

    build_basilisk
else
    print_cyan "Using existing basilisk installation..."
    if [[ "$MODE" == "4" ]]; then
        local_lock_ref="$(read_lock_ref "$BASILISK_DIR/.comphy-lock")"
        local_lock_os="$(read_lock_os "$BASILISK_DIR/.comphy-lock")"
        current_os="$([[ "$OSTYPE" == "darwin"* ]] && echo "darwin" || echo "linux")"

        if [[ -z "$local_lock_ref" ]]; then
            print_red "Error: Existing $BASILISK_DIR is not ref-locked (missing .comphy-lock)."
            print_red "Run with --hard (or remove $BASILISK_DIR) to install ref-locked."
            exit 1
        fi
        if [[ "$local_lock_ref" != "$REF" ]]; then
            print_red "Error: Existing $BASILISK_DIR is locked to ref '$local_lock_ref' but you requested '$REF'."
            print_red "Run with --hard to reinstall the requested ref."
            exit 1
        fi

        # Auto-reinstall on OS mismatch
        if [[ -n "$local_lock_os" ]] && [[ "$local_lock_os" != "$current_os" ]]; then
            print_yellow "OS mismatch detected: installation was built on '$local_lock_os', current OS is '$current_os'."
            print_cyan "Auto-reinstalling for current OS..."
            rm -rf "$BASILISK_DIR"
            install_basilisk_ref_locked "$REF" "$BASILISK_DIR/.comphy-lock"
            build_basilisk
        fi
    fi
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
if [[ "$MODE" == "4" ]]; then
    echo "Lock stamp:"
    echo "  $BASILISK_DIR/.comphy-lock"
fi
