#!/bin/bash
# Linux/macOS version - uses git instead of darcs
# Clones basilisk-source from comphy-lab/basilisk-C GitHub repo
# Ensures that we are always using the latest version from our fork

# Check if --hard flag is passed
HARD_RESET=false
if [[ "$1" == "--hard" ]]; then
    HARD_RESET=true
fi

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

    # Check for git
    if ! command -v git > /dev/null 2>&1; then
        missing_tools+=("git")
    else
        found_tools+=("git")
        printf "\033[0;32m✓ git is installed\033[0m\n"
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
            echo "  brew install gawk"
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

# Function to apply comphy-lab patches (macOS only)
apply_patches() {
    local target_dir="$1"

    if [[ "$OSTYPE" != "darwin"* ]]; then
        # Patches are macOS-specific, skip on other platforms
        return 0
    fi

    printf "\033[0;36mApplying comphy-lab patches...\033[0m\n"

    # Create temp directory for patches
    mkdir -p "$target_dir/.patches_temp"

    # GitHub URLs for patches
    local PATCHES_URL="https://api.github.com/repos/comphy-lab/basilisk-C/contents/patches"
    local RAW_BASE_URL="https://raw.githubusercontent.com/comphy-lab/basilisk-C/main/patches"

    # Get list of patch files (sorted by name for chronological order due to YYYY-MM-DD format)
    local PATCH_FILES
    PATCH_FILES=$(curl -s "$PATCHES_URL" | grep -o '"name": "[^"]*\.patch"' | sed 's/"name": "//;s/"//' | sort)

    if [[ -z "$PATCH_FILES" ]]; then
        printf "\033[0;33mWarning: No patches found or unable to fetch patch list\033[0m\n"
    else
        # Download and apply each patch
        echo "$PATCH_FILES" | while read -r patch_file; do
            if [[ -n "$patch_file" ]]; then
                echo "  Downloading $patch_file..."
                if curl -s -f "$RAW_BASE_URL/$patch_file" -o "$target_dir/.patches_temp/$patch_file"; then
                    echo "  Applying $patch_file..."
                    if (cd "$target_dir" && patch -p1 < ".patches_temp/$patch_file"); then
                        printf "  \033[0;32m✓ Successfully applied $patch_file\033[0m\n"
                    else
                        printf "  \033[0;31m✗ Failed to apply $patch_file\033[0m\n"
                    fi
                else
                    printf "  \033[0;31m✗ Failed to download $patch_file\033[0m\n"
                fi
            fi
        done
    fi

    # Clean up
    rm -rf "$target_dir/.patches_temp"
    echo ""
}

# Function to install basilisk using git sparse checkout from GitHub
install_basilisk() {
    local REPO_URL="https://github.com/comphy-lab/basilisk-C.git"

    printf "\033[0;36mCloning basilisk-source from comphy-lab/basilisk-C (sparse checkout)...\033[0m\n"

    # Use sparse checkout to only get basilisk-source directory
    git clone --depth 1 --filter=blob:none --sparse "$REPO_URL" basilisk-C-temp

    if [[ $? -ne 0 ]]; then
        printf "\033[0;31mError: Failed to clone repository\033[0m\n"
        echo "URL: $REPO_URL"
        exit 1
    fi

    cd basilisk-C-temp || { printf "\033[0;31mError: Failed to change directory to basilisk-C-temp\033[0m\n" >&2; exit 1; }
    git sparse-checkout set basilisk-source
    cd ..

    # Move basilisk-source to basilisk
    printf "\033[0;36mSetting up basilisk directory...\033[0m\n"
    mv basilisk-C-temp/basilisk-source basilisk

    # Clean up temp clone
    rm -rf basilisk-C-temp

    # Apply comphy-lab patches (macOS only)
    apply_patches "basilisk"

    cd basilisk/src || { printf "\033[0;31mError: Failed to change directory to basilisk/src\033[0m\n" >&2; exit 1; }

    if [[ "$OSTYPE" == "darwin"* ]]; then
        printf "\033[0;36mUsing macOS configuration...\033[0m\n"
        ln -s config.osx config
    else
        printf "\033[0;36mUsing Linux configuration...\033[0m\n"
        ln -s config.gcc config
    fi

    printf "\033[0;36mBuilding basilisk (first pass with -k to continue on errors)...\033[0m\n"
    make -k

    printf "\033[0;36mBuilding basilisk (final build)...\033[0m\n"
    make
}

# Check prerequisites first
check_prerequisites

# Remove project config always
rm -rf .project_config

# Check if basilisk needs to be installed
if [[ "$HARD_RESET" == true ]] || [[ ! -d "basilisk" ]]; then
    printf "\033[0;36mInstalling basilisk...\033[0m\n"
    rm -rf basilisk
    install_basilisk
else
    printf "\033[0;36mUsing existing basilisk installation...\033[0m\n"
    cd basilisk/src || { printf "\033[0;31mError: Failed to change directory to basilisk/src\033[0m\n" >&2; exit 1; }
fi

# Setup environment variables
echo "export BASILISK=$PWD" >> ../../.project_config
echo "export PATH=\$PATH:\$BASILISK" >> ../../.project_config

source ../../.project_config

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
