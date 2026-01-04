#!/bin/bash
# Linux version - uses wget and tar instead of darcs (which is broken on Linux)
# Based on https://basilisk.fr/src/INSTALL
# Ensures that we are always using the latest version of basilisk

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
        echo "\033[0;32m✓ make is installed\033[0m"
    fi

    # Check for gawk
    if ! command -v gawk > /dev/null 2>&1; then
        missing_tools+=("gawk")
    else
        found_tools+=("gawk")
        echo "\033[0;32m✓ gawk is installed\033[0m"
    fi

    # Check for wget
    if ! command -v wget > /dev/null 2>&1; then
        missing_tools+=("wget")
    else
        found_tools+=("wget")
        echo "\033[0;32m✓ wget is installed\033[0m"
    fi

    # Check for tar
    if ! command -v tar > /dev/null 2>&1; then
        missing_tools+=("tar")
    else
        found_tools+=("tar")
        echo "\033[0;32m✓ tar is installed\033[0m"
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

        # Linux installation instructions
        echo "To install missing tools on Linux:"
        echo "  sudo apt install ${missing_tools[*]}"
        echo ""

        echo "Please install the missing tools and try again."
        exit 1
    else
        echo "\033[0;32m✅ All prerequisites are satisfied!\033[0m"
        echo ""
    fi
}

# Function to install basilisk using wget
install_basilisk() {
    echo "\033[0;36mDownloading basilisk using wget...\033[0m"
    wget https://basilisk.fr/basilisk/basilisk.tar.gz

    if [[ $? -ne 0 ]]; then
        echo "\033[0;31mError: Failed to download basilisk.tar.gz\033[0m"
        exit 1
    fi

    echo "\033[0;36mExtracting basilisk.tar.gz...\033[0m"
    tar xzf basilisk.tar.gz

    if [[ $? -ne 0 ]]; then
        echo "\033[0;31mError: Failed to extract basilisk.tar.gz\033[0m"
        exit 1
    fi

    # Clean up the tar file
    rm basilisk.tar.gz

    cd basilisk/src

    echo "\033[0;36mUsing Linux configuration...\033[0m"
    ln -s config.gcc config

    echo "\033[0;36mBuilding basilisk (first pass with -k to continue on errors)...\033[0m"
    make -k

    echo "\033[0;36mBuilding basilisk (final build)...\033[0m"
    make
}

# Check prerequisites first
check_prerequisites

# Remove project config always
rm -rf .project_config

# Check if basilisk needs to be installed
if [[ "$HARD_RESET" == true ]] || [[ ! -d "basilisk" ]]; then
    echo "\033[0;36mInstalling basilisk...\033[0m"
    rm -rf basilisk
    install_basilisk
else
    echo "\033[0;36mUsing existing basilisk installation...\033[0m"
    cd basilisk/src
fi

# Setup environment variables
echo "export BASILISK=$PWD" >> ../../.project_config
echo "export PATH=\$PATH:\$BASILISK" >> ../../.project_config

source ../../.project_config

# Check if qcc is working properly
echo ""
echo "\033[0;36mChecking qcc installation...\033[0m"
if ! qcc --version > /dev/null 2>&1; then
    echo "\033[0;31mError: qcc is not working properly.\033[0m"
    echo "Please ensure you have build-essential installed."
    echo "You can install it by running: sudo apt install build-essential"
    echo "For more details, visit: https://basilisk.fr/src/INSTALL"
    exit 1
else
    echo "\033[0;32m✅ qcc is properly installed.\033[0m"
    qcc --version
fi