#!/bin/zsh
# tested on MacOS only. Let us know if you find issues running with Linux by opening an issue. 
# modify using https://basilisk.fr/src/INSTALL 
# ensures that we are always using the latest version of basilisk

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

# Function to install basilisk
install_basilisk() {
    if ! darcs clone https://basilisk.fr/basilisk "$BASILISK_DIR"; then
        echo "\033[0;31mError: Failed to clone basilisk into $BASILISK_DIR\033[0m"
        exit 1
    fi
    if ! cd "$BASILISK_DIR"; then
        echo "\033[0;31mError: Failed to change directory to $BASILISK_DIR\033[0m"
        exit 1
    fi

    # Apply comphy-lab patches (macOS only)
    if [[ "$OSTYPE" == "darwin"* ]]; then
        echo "Downloading and applying comphy-lab patches..."
        local patch_failed=false

        # Create temp directory for patches
        mkdir -p .patches_temp

        # GitHub URLs for patches
        PATCHES_URL="https://api.github.com/repos/comphy-lab/basilisk-C/contents/patches"
        RAW_BASE_URL="https://raw.githubusercontent.com/comphy-lab/basilisk-C/main/patches"

        # Get list of patch files (sorted by name which gives chronological order due to YYYY-MM-DD format)
        PATCH_FILES=$(curl -s "$PATCHES_URL" | grep -o '"name": "[^"]*\.patch"' | sed 's/"name": "//;s/"//' | sort)

        if [[ -z "$PATCH_FILES" ]]; then
            echo "\033[0;33mWarning: No patches found or unable to fetch patch list\033[0m"
        else
            # Download and apply each patch
            while read -r patch_file; do
                if [[ -n "$patch_file" ]]; then
                    # Skip local-bview patch unless --local-bview flag was provided
                    if [[ "$patch_file" == *"-local-bview.patch" ]] && [[ "$LOCAL_BVIEW" != "true" ]]; then
                        echo "  Skipping $patch_file (use --local-bview to apply)"
                        continue
                    fi
                    echo "  Downloading $patch_file..."
                    if curl -s -f "$RAW_BASE_URL/$patch_file" -o ".patches_temp/$patch_file"; then
                        echo "  Applying $patch_file..."
                        if patch -p1 < ".patches_temp/$patch_file"; then
                            echo "  \033[0;32m✓ Successfully applied $patch_file\033[0m"
                        else
                            echo "  \033[0;31m✗ Failed to apply $patch_file\033[0m"
                            patch_failed=true
                        fi
                    else
                        echo "  \033[0;31m✗ Failed to download $patch_file\033[0m"
                        patch_failed=true
                    fi
                fi
            done <<< "$PATCH_FILES"
        fi

        # Clean up
        rm -rf .patches_temp
        echo ""

        if [[ "$patch_failed" == true ]]; then
            echo "\033[0;31mError: Failed to apply comphy-lab patches in $BASILISK_DIR\033[0m"
            exit 1
        fi
    fi

    if ! cd "$BASILISK_SRC_DIR"; then
        echo "\033[0;31mError: Failed to change directory to $BASILISK_SRC_DIR\033[0m"
        exit 1
    fi

    if [[ "$OSTYPE" == "darwin"* ]]; then
        echo "Using MacOS"
        ln -s config.osx config
    else
        echo "Using Linux"
        ln -s config.gcc config
    fi
    if ! make -k; then
        echo "\033[0;31mError: make -k failed in $BASILISK_SRC_DIR\033[0m"
        exit 1
    fi
    if ! make; then
        echo "\033[0;31mError: make failed in $BASILISK_SRC_DIR\033[0m"
        exit 1
    fi
}

# Check prerequisites first
check_prerequisites

# Remove project config always
rm -rf "$PROJECT_CONFIG"

# Check if basilisk needs to be installed
if [[ "$HARD_RESET" == true ]] || [[ ! -d "$BASILISK_DIR" ]]; then
    echo "Installing basilisk..."
    rm -rf "$BASILISK_DIR"
    install_basilisk
else
    echo "Using existing basilisk installation..."
    if [[ ! -d "$BASILISK_SRC_DIR" ]]; then
        echo "\033[0;31mError: Missing $BASILISK_SRC_DIR. Run with --hard to reinstall.\033[0m"
        exit 1
    fi
    if ! cd "$BASILISK_SRC_DIR"; then
        echo "\033[0;31mError: Failed to change directory to $BASILISK_SRC_DIR\033[0m"
        exit 1
    fi
fi

# Setup environment variables
if [[ ! -d "$BASILISK_SRC_DIR" ]]; then
    echo "\033[0;31mError: Expected $BASILISK_SRC_DIR to exist before writing .project_config\033[0m"
    exit 1
fi
printf "export BASILISK=%s\n" "$BASILISK_SRC_DIR" > "$PROJECT_CONFIG"
# Prepend BASILISK to PATH so Basilisk tools (qcc, etc.) take precedence
printf "export PATH=\\$BASILISK:\\$PATH\n" >> "$PROJECT_CONFIG"

source "$PROJECT_CONFIG"

# Check if qcc is working properly
printf '\nChecking qcc installation...\n'
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
