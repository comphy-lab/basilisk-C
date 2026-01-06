#!/bin/zsh
# Local patch testing script - applies patches from local patches/ directory
# Use this to test patches before pushing to GitHub

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
PATCHES_DIR="$SCRIPT_DIR/patches"

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
        echo "\033[0;32m✓ make is installed\033[0m"
    fi

    # Check for darcs
    if ! command -v darcs > /dev/null 2>&1; then
        missing_tools+=("darcs")
    else
        found_tools+=("darcs")
        echo "\033[0;32m✓ darcs is installed\033[0m"
    fi

    # Check for gcc
    if ! command -v gcc > /dev/null 2>&1; then
        missing_tools+=("gcc")
    else
        found_tools+=("gcc")
        echo "\033[0;32m✓ gcc is installed\033[0m"
    fi

    echo ""

    if [[ ${#missing_tools[@]} -gt 0 ]]; then
        echo "\033[0;31mError: Missing required tools: ${missing_tools[*]}\033[0m"
        echo ""

        if [[ "$OSTYPE" == "darwin"* ]]; then
            # macOS installation instructions
            for tool in "${missing_tools[@]}"; do
                case "$tool" in
                    "make")
                        echo "To install make:"
                        echo "  xcode-select --install"
                        echo ""
                        ;;
                    "darcs")
                        echo "To install darcs:"
                        echo "  brew install darcs"
                        echo ""
                        ;;
                    "gcc")
                        echo "To install gcc:"
                        echo "  xcode-select --install"
                        echo ""
                        ;;
                esac
            done
        else
            # Linux installation instructions
            for tool in "${missing_tools[@]}"; do
                case "$tool" in
                    "make")
                        echo "To install make on Linux:"
                        echo "  Ubuntu/Debian: sudo apt-get install build-essential"
                        echo "  RHEL/CentOS: sudo yum groupinstall 'Development Tools'"
                        echo ""
                        ;;
                    "darcs")
                        echo "To install darcs:"
                        echo "  Visit https://darcs.net/ for installation instructions"
                        echo ""
                        ;;
                    "gcc")
                        echo "To install gcc on Linux:"
                        echo "  Ubuntu/Debian: sudo apt-get install build-essential"
                        echo "  RHEL/CentOS: sudo yum groupinstall 'Development Tools'"
                        echo ""
                        ;;
                esac
            done
        fi

        echo "Please install the missing tools and try again."
        exit 1
    else
        echo "\033[0;32m✅ All prerequisites are satisfied!\033[0m"
        echo ""
    fi
}

# Function to apply patches from local patches/ directory
apply_local_patches() {
    local target_dir="$1"
    local apply_local_bview="$2"

    if [[ ! -d "$PATCHES_DIR" ]]; then
        echo "\033[0;31mError: Patches directory not found: $PATCHES_DIR\033[0m"
        return 1
    fi

    printf "\033[0;36mApplying patches from local directory: $PATCHES_DIR\033[0m\n"

    # Get list of patch files sorted by name (YYYY-MM-DD format ensures chronological order)
    local patch_files=($(ls "$PATCHES_DIR"/*.patch 2>/dev/null | sort))

    if [[ ${#patch_files[@]} -eq 0 ]]; then
        printf "\033[0;33mWarning: No patches found in $PATCHES_DIR\033[0m\n"
        return 0
    fi

    for patch_file in "${patch_files[@]}"; do
        local patch_name=$(basename "$patch_file")

        # Skip local-bview patch unless --local-bview flag was provided
        if [[ "$patch_name" == *"local-bview"* ]] && [[ "$apply_local_bview" != "true" ]]; then
            echo "  Skipping $patch_name (use --local-bview to apply)"
            continue
        fi

        echo "  Applying $patch_name..."
        if (cd "$target_dir" && patch -p1 < "$patch_file"); then
            printf "  \033[0;32m✓ Successfully applied $patch_name\033[0m\n"
        else
            printf "  \033[0;31m✗ Failed to apply $patch_name\033[0m\n"
            printf "  \033[0;33mCheck $target_dir/src/*.rej for details\033[0m\n"
        fi
    done
    echo ""
}

# Function to install basilisk
install_basilisk() {
    darcs clone https://basilisk.fr/basilisk
    cd basilisk

    # Apply patches from local patches/ directory
    apply_local_patches "$(pwd)" "$LOCAL_BVIEW"

    cd src

    if [[ "$OSTYPE" == "darwin"* ]]; then
        echo "Using MacOS"
        ln -s config.osx config
    else
        echo "Using Linux"
        ln -s config.gcc config
    fi
    make -k
    make
}

# Check prerequisites first
check_prerequisites

# Remove project config always
rm -rf .project_config

# Check if basilisk needs to be installed
if [[ "$HARD_RESET" == true ]] || [[ ! -d "basilisk" ]]; then
    echo "Installing basilisk..."
    rm -rf basilisk
    install_basilisk
else
    echo "Using existing basilisk installation..."
    cd basilisk/src
fi

# Setup environment variables
echo "export BASILISK=$PWD" >> ../../.project_config
echo "export PATH=\$BASILISK:\$PATH" >> ../../.project_config

source ../../.project_config

# Check if qcc is working properly
echo "\nChecking qcc installation..."
if ! qcc --version > /dev/null 2>&1; then
    echo "\033[0;31mError: qcc is not working properly.\033[0m"
    echo "Please ensure you have Xcode Command Line Tools installed."
    echo "You can install them by running: xcode-select --install"
    echo "For more details, visit: http://basilisk.fr/src/INSTALL"
    exit 1
else
    echo "\033[0;32mqcc is properly installed.\033[0m"
    qcc --version
fi
