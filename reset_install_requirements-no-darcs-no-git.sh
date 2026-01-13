#!/bin/bash
# Linux/macOS version - uses wget and tar instead of darcs
# Based on https://basilisk.fr/src/INSTALL
# Ensures that we are always using the latest version of basilisk from basilisk.fr

# Check for flags
HARD_RESET=false
LOCAL_BVIEW=false

for arg in "$@"; do
    case "$arg" in
        --hard)
            HARD_RESET=true
            ;;
        --local-bview)
            LOCAL_BVIEW=true
            ;;
    esac
done

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$SCRIPT_DIR"
PROJECT_CONFIG="$REPO_ROOT/.project_config"
BASILISK_DIR="$REPO_ROOT/basilisk"
BASILISK_SRC_DIR="$BASILISK_DIR/src"
PATCHES_DIR="$REPO_ROOT/patches"

# Function to check prerequisites
check_prerequisites() {
    local missing_tools=()
    local found_tools=()

    echo "Checking prerequisites..."
    echo ""

    # Check for make
    if ! command -v make > /dev/null 2>&1; then
        missing_tools+=("make")
    else
        found_tools+=("make")
        printf "\033[0;32m✓ make is installed\033[0m\n"
    fi

    # Check for gawk
    if ! command -v gawk > /dev/null 2>&1; then
        missing_tools+=("gawk")
    else
        found_tools+=("gawk")
        printf "\033[0;32m✓ gawk is installed\033[0m\n"
    fi

    # Check for wget
    if ! command -v wget > /dev/null 2>&1; then
        missing_tools+=("wget")
    else
        found_tools+=("wget")
        printf "\033[0;32m✓ wget is installed\033[0m\n"
    fi

    # Check for curl (needed for applying patches on macOS)
    if ! command -v curl > /dev/null 2>&1; then
        if [[ "$OSTYPE" == "darwin"* ]]; then
            missing_tools+=("curl")
        fi
    else
        found_tools+=("curl")
        printf "\033[0;32m✓ curl is installed\033[0m\n"
    fi

    # Check for tar
    if ! command -v tar > /dev/null 2>&1; then
        missing_tools+=("tar")
    else
        found_tools+=("tar")
        printf "\033[0;32m✓ tar is installed\033[0m\n"
    fi

    # Check for gcc
    if ! command -v gcc > /dev/null 2>&1; then
        missing_tools+=("gcc")
    else
        found_tools+=("gcc")
        printf "\033[0;32m✓ gcc is installed\033[0m\n"
    fi

    echo ""

    if [[ ${#missing_tools[@]} -gt 0 ]]; then
        printf "\033[0;31mError: Missing required tools: ${missing_tools[*]}\033[0m\n"
        echo ""

        if [[ "$OSTYPE" == "darwin"* ]]; then
            # macOS installation instructions
            echo "To install missing tools on macOS:"
            echo "  xcode-select --install"
            echo "  brew install gawk wget"
        else
            # Linux installation instructions
            echo "To install missing tools on Linux:"
            echo "  sudo apt install ${missing_tools[*]}"
        fi
        echo ""

        echo "Please install the missing tools and try again."
        exit 1
    else
        printf "\033[0;32m✅ All prerequisites are satisfied!\033[0m\n"
        echo ""
    fi
}

# Function to apply comphy-lab patches
# Checks for local patches/ directory first, falls back to GitHub if not available
apply_patches() {
    local target_dir="$1"
    local apply_local_bview="${2:-false}"
    local patch_failed=false

    printf "\033[0;36mApplying comphy-lab patches...\033[0m\n"

    # Check for local patches directory first (useful for HPC/offline environments)
    if [[ -d "$PATCHES_DIR" ]]; then
        printf "\033[0;36mUsing local patches from: $PATCHES_DIR\033[0m\n"

        # Get list of patch files sorted by name (YYYY-MM-DD format ensures chronological order)
        local patch_files
        patch_files=($(ls "$PATCHES_DIR"/*.patch 2>/dev/null | sort))

        if [[ ${#patch_files[@]} -eq 0 ]]; then
            printf "\033[0;33mWarning: No patches found in $PATCHES_DIR\033[0m\n"
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
                    printf "  \033[0;32m✓ Successfully applied $patch_name\033[0m\n"
                else
                    printf "  \033[0;31m✗ Failed to apply $patch_name\033[0m\n"
                    patch_failed=true
                fi
            done
        fi
    else
        # Fall back to downloading patches from GitHub (requires internet)
        printf "\033[0;36mNo local patches/ directory found, downloading from GitHub...\033[0m\n"
        printf "\033[0;33mNote: Internet connection required. For offline use, include patches/ directory.\033[0m\n"

        # Create temp directory for patches
        mkdir -p "$target_dir/.patches_temp"

        # GitHub URLs for patches
        local PATCHES_URL="https://api.github.com/repos/comphy-lab/basilisk-C/contents/patches"
        local RAW_BASE_URL="https://raw.githubusercontent.com/comphy-lab/basilisk-C/main/patches"

        # Get list of patch files (sorted by name for chronological order due to YYYY-MM-DD format)
        local PATCH_FILES
        PATCH_FILES=$(curl -s "$PATCHES_URL" | grep -o '"name": "[^"]*\.patch"' | sed 's/"name": "//;s/"//' | sort)

        if [[ -z "$PATCH_FILES" ]]; then
            printf "\033[0;33mWarning: No patches found or unable to fetch patch list (check internet connection)\033[0m\n"
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
                    if curl -s -f "$RAW_BASE_URL/$patch_file" -o "$target_dir/.patches_temp/$patch_file"; then
                        echo "  Applying $patch_file..."
                        if (cd "$target_dir" && patch -p1 < ".patches_temp/$patch_file"); then
                            printf "  \033[0;32m✓ Successfully applied $patch_file\033[0m\n"
                        else
                            printf "  \033[0;31m✗ Failed to apply $patch_file\033[0m\n"
                            patch_failed=true
                        fi
                    else
                        printf "  \033[0;31m✗ Failed to download $patch_file\033[0m\n"
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

# Function to install basilisk using wget
install_basilisk() {
    if ! cd "$REPO_ROOT"; then
        printf "\033[0;31mError: Failed to change directory to $REPO_ROOT\033[0m\n" >&2
        exit 1
    fi
    printf "\033[0;36mDownloading basilisk using wget...\033[0m\n"
    wget https://basilisk.fr/basilisk/basilisk.tar.gz

    if [[ $? -ne 0 ]]; then
        printf "\033[0;31mError: Failed to download basilisk.tar.gz\033[0m\n"
        exit 1
    fi

    printf "\033[0;36mExtracting basilisk.tar.gz...\033[0m\n"
    tar xzf basilisk.tar.gz

    if [[ $? -ne 0 ]]; then
        printf "\033[0;31mError: Failed to extract basilisk.tar.gz\033[0m\n"
        exit 1
    fi

    # Clean up the tar file
    rm basilisk.tar.gz

    # Apply comphy-lab patches (macOS only)
    if ! apply_patches "$BASILISK_DIR" "$LOCAL_BVIEW"; then
        printf "\033[0;31mError: Failed to apply comphy-lab patches in $BASILISK_DIR\033[0m\n" >&2
        exit 1
    fi

    if ! cd "$BASILISK_SRC_DIR"; then
        printf "\033[0;31mError: Failed to change directory to $BASILISK_SRC_DIR\033[0m\n" >&2
        exit 1
    fi

    if [[ "$OSTYPE" == "darwin"* ]]; then
        printf "\033[0;36mUsing macOS configuration...\033[0m\n"
        ln -s config.osx config
    else
        printf "\033[0;36mUsing Linux configuration...\033[0m\n"
        ln -s config.gcc config
    fi

    printf "\033[0;36mBuilding basilisk (first pass with -k to continue on errors)...\033[0m\n"
    if ! make -k; then
        printf "\033[0;31mError: make -k failed in $BASILISK_SRC_DIR\033[0m\n" >&2
        exit 1
    fi

    printf "\033[0;36mBuilding basilisk (final build)...\033[0m\n"
    if ! make; then
        printf "\033[0;31mError: make failed in $BASILISK_SRC_DIR\033[0m\n" >&2
        exit 1
    fi
}

# Check prerequisites first
check_prerequisites

# Remove project config always
rm -rf "$PROJECT_CONFIG"

# Check if basilisk needs to be installed
if [[ "$HARD_RESET" == true ]] || [[ ! -d "$BASILISK_DIR" ]]; then
    printf "\033[0;36mInstalling basilisk...\033[0m\n"
    rm -rf "$BASILISK_DIR"
    install_basilisk
else
    printf "\033[0;36mUsing existing basilisk installation...\033[0m\n"
    if [[ ! -d "$BASILISK_SRC_DIR" ]]; then
        printf "\033[0;31mError: Missing $BASILISK_SRC_DIR. Run with --hard to reinstall.\033[0m\n" >&2
        exit 1
    fi
    if ! cd "$BASILISK_SRC_DIR"; then
        printf "\033[0;31mError: Failed to change directory to $BASILISK_SRC_DIR\033[0m\n" >&2
        exit 1
    fi
fi

# Setup environment variables
if [[ ! -d "$BASILISK_SRC_DIR" ]]; then
    printf "\033[0;31mError: Expected $BASILISK_SRC_DIR to exist before writing .project_config\033[0m\n" >&2
    exit 1
fi
printf "export BASILISK=%s\n" "$BASILISK_SRC_DIR" > "$PROJECT_CONFIG"
# Prepend BASILISK to PATH so Basilisk tools (qcc, etc.) take precedence
printf "export PATH=\$BASILISK:\$PATH\n" >> "$PROJECT_CONFIG"

source "$PROJECT_CONFIG"

# Check if qcc is working properly
echo ""
printf "\033[0;36mChecking qcc installation...\033[0m\n"
if ! qcc --version > /dev/null 2>&1; then
    printf "\033[0;31mError: qcc is not working properly.\033[0m\n"
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
    printf "\033[0;32m✅ qcc is properly installed.\033[0m\n"
    qcc --version
fi
