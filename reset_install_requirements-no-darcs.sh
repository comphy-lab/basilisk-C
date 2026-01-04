#!/bin/bash
# Linux/macOS version - uses wget and tar instead of darcs
# Downloads basilisk-source from comphy-lab/basilisk-C GitHub repo
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
        echo "\033[0;32m✓ make is installed\033[0m"
    fi

    # Check for gawk
    if ! command -v gawk > /dev/null 2>&1; then
        missing_tools+=("gawk")
    else
        found_tools+=("gawk")
        echo "\033[0;32m✓ gawk is installed\033[0m"
    fi

    # Check for wget or curl
    if command -v wget > /dev/null 2>&1; then
        found_tools+=("wget")
        echo "\033[0;32m✓ wget is installed\033[0m"
        USE_WGET=true
    elif command -v curl > /dev/null 2>&1; then
        found_tools+=("curl")
        echo "\033[0;32m✓ curl is installed\033[0m"
        USE_WGET=false
    else
        missing_tools+=("wget or curl")
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
        echo "\033[0;32m✅ All prerequisites are satisfied!\033[0m"
        echo ""
    fi
}

# Function to download file using wget or curl
download_file() {
    local url="$1"
    local output="$2"

    if [[ "$USE_WGET" == true ]]; then
        wget -O "$output" "$url"
    else
        curl -L -o "$output" "$url"
    fi
}

# Function to install basilisk using wget/curl from GitHub
install_basilisk() {
    local DOWNLOAD_URL="https://raw.githubusercontent.com/comphy-lab/basilisk-C/main/basilisk-source.tar.gz"

    echo "\033[0;36mDownloading basilisk-source from comphy-lab/basilisk-C...\033[0m"
    download_file "$DOWNLOAD_URL" "basilisk-source.tar.gz"

    if [[ $? -ne 0 ]]; then
        echo "\033[0;31mError: Failed to download basilisk-source.tar.gz\033[0m"
        echo "URL: $DOWNLOAD_URL"
        exit 1
    fi

    echo "\033[0;36mExtracting basilisk-source.tar.gz...\033[0m"
    tar xzf basilisk-source.tar.gz

    if [[ $? -ne 0 ]]; then
        echo "\033[0;31mError: Failed to extract basilisk-source.tar.gz\033[0m"
        exit 1
    fi

    # Clean up the tar file
    rm basilisk-source.tar.gz

    cd basilisk-source/src

    if [[ "$OSTYPE" == "darwin"* ]]; then
        echo "\033[0;36mUsing macOS configuration...\033[0m"
        ln -s config.osx config
    else
        echo "\033[0;36mUsing Linux configuration...\033[0m"
        ln -s config.gcc config
    fi

    echo "\033[0;36mBuilding basilisk (first pass with -k to continue on errors)...\033[0m"
    make -k

    echo "\033[0;36mBuilding basilisk (final build)...\033[0m"
    make
}

# Check prerequisites first
check_prerequisites

# Remove project config always
rm -rf .project_config

# Check if basilisk-source needs to be installed
if [[ "$HARD_RESET" == true ]] || [[ ! -d "basilisk-source" ]]; then
    echo "\033[0;36mInstalling basilisk-source...\033[0m"
    rm -rf basilisk-source
    install_basilisk
else
    echo "\033[0;36mUsing existing basilisk-source installation...\033[0m"
    cd basilisk-source/src
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
    echo "\033[0;32m✅ qcc is properly installed.\033[0m"
    qcc --version
fi
