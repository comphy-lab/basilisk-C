#!/bin/bash
# Ref-locked Basilisk installation script
# Installs Basilisk from OS-specific tarballs attached to a GitHub Release tag
# in comphy-lab/basilisk-C for reproducible installs (no darcs/git required).
#
# Supported remote usage:
#   curl -sL https://raw.githubusercontent.com/comphy-lab/basilisk-C/main/reset_install_basilisk-ref-locked.sh | bash -s -- --ref=<tag>
#   curl -sL https://raw.githubusercontent.com/comphy-lab/basilisk-C/main/reset_install_basilisk-ref-locked.sh | zsh  -s -- --ref=<tag>

set -euo pipefail

HARD_RESET=false
LOCAL_BVIEW=false
SHOW_HELP=false
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
    --ref=*)
      REF="${arg#*=}"
      ;;
  esac
done

print_green() { printf "\033[0;32m%s\033[0m\n" "$1"; }
print_red() { printf "\033[0;31m%s\033[0m\n" "$1"; }
print_yellow() { printf "\033[0;33m%s\033[0m\n" "$1"; }
print_cyan() { printf "\033[0;36m%s\033[0m\n" "$1"; }

PATCHES_API_URL="https://api.github.com/repos/comphy-lab/basilisk-C/contents/patches"

show_help() {
  cat << 'EOF'
Ref-locked Basilisk Installation Script
======================================

Installs Basilisk from OS-specific tarballs attached to a GitHub Release tag in
comphy-lab/basilisk-C.

Usage:
  ./reset_install_basilisk-ref-locked.sh [OPTIONS]

Options:
  --help, -h      Show this help message
  --hard          Force reinstall (removes existing basilisk directory)
  --local-bview   Apply optional local-bview convenience patch (enables `bview --local` URL output for bview-local-client)
  --ref=REF       Required. GitHub Release tag in comphy-lab/basilisk-C

Examples:
  ./reset_install_basilisk-ref-locked.sh --ref=v2026-01-13 --hard
  ./reset_install_basilisk-ref-locked.sh --ref=v2026-01-13 --local-bview --hard

Remote:
  curl -sL https://raw.githubusercontent.com/comphy-lab/basilisk-C/main/reset_install_basilisk-ref-locked.sh | bash -s -- --ref=v2026-01-13 --hard
  curl -sL https://raw.githubusercontent.com/comphy-lab/basilisk-C/main/reset_install_basilisk-ref-locked.sh | zsh  -s -- --ref=v2026-01-13 --hard

Notes:
  - GitHub Release tarballs intentionally exclude the local-bview patch; `--local-bview` downloads and applies it for the same `--ref`.
EOF
}

check_tool() {
  local tool="$1"
  if command -v "$tool" >/dev/null 2>&1; then
    print_green "✓ $tool is installed"
    return 0
  fi
  return 1
}

show_install_instructions() {
  local tool="$1"
  if [[ "$OSTYPE" == "darwin"* ]]; then
    case "$tool" in
      make|gcc)
        echo "  xcode-select --install"
        ;;
      gawk)
        echo "  brew install gawk"
        ;;
      curl)
        echo "  brew install curl"
        ;;
      tar|patch)
        echo "  (should be pre-installed on macOS)"
        ;;
    esac
  else
    case "$tool" in
      make|gcc)
        echo "  sudo apt install build-essential"
        ;;
      gawk|curl|tar|patch)
        echo "  sudo apt install $tool"
        ;;
    esac
  fi
}

check_prerequisites() {
  local missing_tools=()

  echo "Checking prerequisites..."
  echo ""

  check_tool "make" || missing_tools+=("make")
  check_tool "gcc" || missing_tools+=("gcc")
  check_tool "curl" || missing_tools+=("curl")
  check_tool "tar" || missing_tools+=("tar")
  check_tool "gawk" || missing_tools+=("gawk")
  if [[ "$LOCAL_BVIEW" == "true" ]]; then
    check_tool "patch" || missing_tools+=("patch")
  fi

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
  fi

  print_green "✅ All prerequisites are satisfied!"
  echo ""
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
    printf "created_utc=%s\n" "$(date -u +"%Y-%m-%dT%H:%M:%SZ" 2>/dev/null || date)"
    printf "local_bview=%s\n" "$apply_local_bview"
    printf "patches_applied=%s\n" "${applied_patches[*]}"
    printf "patches_skipped=%s\n" "${skipped_patches[*]}"
  } > "$lock_file"
}

build_basilisk() {
  if ! cd "$BASILISK_SRC_DIR"; then
    print_red "Error: Failed to change directory to $BASILISK_SRC_DIR"
    exit 1
  fi

  if [[ -e config ]]; then
    rm -f config
  fi

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

write_project_config() {
  local project_config="$1"
  local basilisk_src_dir="$2"

  {
    printf 'export BASILISK="%s"\n' "$basilisk_src_dir"
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
  } > "$project_config"
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

  if ! qcc --version >/dev/null 2>&1; then
    print_red "Error: qcc is not working properly."
    if [[ "$OSTYPE" == "darwin"* ]]; then
      echo "Please ensure you have Xcode Command Line Tools installed."
      echo "You can install them by running: xcode-select --install"
    else
      echo "Please ensure you have build-essential installed."
      echo "You can install it by running: sudo apt install build-essential"
    fi
    exit 1
  fi

  print_green "✅ qcc is properly installed."
  qcc --version
}

if [[ "$SHOW_HELP" == true ]]; then
  show_help
  exit 0
fi

if [[ -z "${REF}" ]]; then
  print_red "Error: --ref is required"
  echo ""
  show_help
  exit 1
fi

REPO_ROOT="$PWD"
if [[ -n "${0:-}" ]] && [[ -f "${0:-}" ]]; then
  SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
  REPO_ROOT="$SCRIPT_DIR"
fi

PROJECT_CONFIG="$REPO_ROOT/.project_config"
BASILISK_DIR="$REPO_ROOT/basilisk"
BASILISK_SRC_DIR="$BASILISK_DIR/src"
PATCHES_DIR="$REPO_ROOT/patches"
LOCK_FILE="$BASILISK_DIR/.comphy-lock"

check_prerequisites

rm -rf "$PROJECT_CONFIG"

if [[ "$HARD_RESET" == true ]] || [[ ! -d "$BASILISK_DIR" ]]; then
  print_cyan "Installing basilisk (ref-locked)..."
  rm -rf "$BASILISK_DIR"

  if [[ "$OSTYPE" == "darwin"* ]]; then
    asset_name="basilisk-mac.tar.gz"
  else
    asset_name="basilisk-linux.tar.gz"
  fi
  download_url="https://github.com/comphy-lab/basilisk-C/releases/download/${REF}/${asset_name}"

  temp_dir="$(mktemp -d 2>/dev/null || mktemp -d -t basilisk-C)"
  if [[ -z "$temp_dir" ]] || [[ ! -d "$temp_dir" ]]; then
    print_red "Error: Failed to create temp directory"
    exit 1
  fi

  print_cyan "Downloading: $download_url"
  if ! curl -s -f -L "$download_url" -o "$temp_dir/$asset_name"; then
    print_red "Error: Failed to download release asset for ref '$REF'"
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

  patches_dest=""
  if [[ -d "$temp_dir/patches" ]]; then
    if [[ ! -d "$PATCHES_DIR" ]]; then
      if ! mv "$temp_dir/patches" "$PATCHES_DIR"; then
        print_red "Error: Failed to move patches into $PATCHES_DIR"
        rm -rf "$temp_dir"
        exit 1
      fi
      patches_dest="$PATCHES_DIR"
    else
      patches_dest="$BASILISK_DIR/.comphy-patches"
      if ! mv "$temp_dir/patches" "$patches_dest"; then
        print_red "Error: Failed to move patches into $patches_dest"
        rm -rf "$temp_dir"
        exit 1
      fi
    fi
  else
    print_red "Error: Tarball is missing expected patches/ directory"
    rm -rf "$temp_dir"
    exit 1
  fi

  if [[ "$LOCAL_BVIEW" == true ]]; then
    print_cyan "Applying local-bview patch (from pinned ref)..."

    patches_api_url_ref="${PATCHES_API_URL}?ref=${REF}"
    raw_base_url="https://raw.githubusercontent.com/comphy-lab/basilisk-C/${REF}/patches"
    patch_file="$(curl -s "$patches_api_url_ref" | grep -o '\"name\": \"[^\"]*\.patch\"' | sed 's/\"name\": \"//;s/\"//' | grep -E -- '-local-bview\.patch$' | head -n 1)"

    if [[ -z "$patch_file" ]]; then
      print_red "Error: Could not find a *-local-bview.patch file at ref '$REF'"
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

  write_lock_stamp "$LOCK_FILE" "$REF" "$patches_dest" "$LOCAL_BVIEW"
  rm -rf "$temp_dir"

  build_basilisk
else
  print_cyan "Using existing basilisk installation..."

  local_lock_ref="$(read_lock_ref "$LOCK_FILE")"
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

  if [[ ! -d "$BASILISK_SRC_DIR" ]]; then
    print_red "Error: Missing $BASILISK_SRC_DIR. Run with --hard to reinstall."
    exit 1
  fi
fi

write_project_config "$PROJECT_CONFIG" "$BASILISK_SRC_DIR"
source "$PROJECT_CONFIG"
verify_installation

echo ""
print_green "✅ Basilisk ref-locked installation complete!"
echo "To use basilisk in your shell, run:"
echo "  source $PROJECT_CONFIG"
echo "Lock stamp:"
echo "  $LOCK_FILE"
