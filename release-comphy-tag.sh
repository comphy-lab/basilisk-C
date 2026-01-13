#!/bin/bash
# Create a comphy-lab Basilisk release tag and publish OS-specific, pre-patched
# Basilisk source tarballs as GitHub Release assets.
#
# - Source snapshot comes from upstream (darcs) pulled into basilisk-source/
# - Patches come from patches/ and are applied before packaging
# - Produces two tarballs:
#     basilisk-mac.tar.gz   (includes macOS patches + all non-macOS patches)
#     basilisk-linux.tar.gz (includes all non-macOS patches)
#   The local-bview patch is NOT applied to either tarball.
#
# This script also runs install tests:
# - On macOS: tests macOS install natively and Linux install in Docker
# - On Linux: tests Linux install natively

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$SCRIPT_DIR"

TAG=""
DRY_RUN=false
SKIP_TESTS=false
FORCE=false

usage() {
  cat <<'EOF'
Usage: ./release-comphy-tag.sh [OPTIONS]

Options:
  --tag=TAG        Release tag name (default: vYYYY-MM-DD in UTC)
  --skip-tests     Skip install tests
  --dry-run        Print actions without pushing/tagging/releasing
  --force          Delete/recreate local tag if it exists (still fails if remote tag exists)
  --help, -h       Show this help message

Notes:
  - Requires gh auth and push access to the repo.
  - Produces and uploads these GitHub Release assets:
      basilisk-mac.tar.gz (+ .sha256)
      basilisk-linux.tar.gz (+ .sha256)
EOF
}

for arg in "$@"; do
  case "$arg" in
    --tag=*)
      TAG="${arg#*=}"
      ;;
    --skip-tests)
      SKIP_TESTS=true
      ;;
    --dry-run)
      DRY_RUN=true
      ;;
    --force)
      FORCE=true
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "Error: Unknown argument: $arg" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if [[ -z "$TAG" ]]; then
  TAG="v$(date -u +%Y-%m-%d)"
fi

print_cyan() { printf "\033[0;36m%s\033[0m\n" "$1"; }
print_green() { printf "\033[0;32m%s\033[0m\n" "$1"; }
print_red() { printf "\033[0;31m%s\033[0m\n" "$1"; }

require_tool() {
  local tool="$1"
  if ! command -v "$tool" >/dev/null 2>&1; then
    print_red "Error: Missing required tool: $tool"
    exit 1
  fi
}

sha256_file() {
  local file="$1"
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$file"
  else
    shasum -a 256 "$file"
  fi
}

copy_basilisk_source_to_stage() {
  local stage_basilisk_dir="$1"
  mkdir -p "$stage_basilisk_dir"
  tar -C "$REPO_ROOT/basilisk-source" -cf - --exclude "_darcs" . | tar -C "$stage_basilisk_dir" -xf -
}

copy_selected_patches() {
  local platform="$1" # mac | linux
  local dest_patches_dir="$2"
  mkdir -p "$dest_patches_dir"

  local patch_file
  for patch_file in "$REPO_ROOT"/patches/*.patch; do
    local patch_name
    patch_name="$(basename "$patch_file")"

    if [[ "$patch_name" == *"-local-bview.patch" ]]; then
      continue
    fi
    if [[ "$platform" == "linux" ]] && [[ "$patch_name" == *"-macos-"* ]]; then
      continue
    fi

    cp "$patch_file" "$dest_patches_dir/"
  done
}

apply_patches_in_dir() {
  local target_dir="$1"
  local patches_dir="$2"

  local patch_files=()
  patch_files=($(ls "$patches_dir"/*.patch 2>/dev/null | sort))
  if [[ ${#patch_files[@]} -eq 0 ]]; then
    print_red "Error: No patches found in $patches_dir"
    exit 1
  fi

  local patch_path
  for patch_path in "${patch_files[@]}"; do
    local patch_name
    patch_name="$(basename "$patch_path")"
    print_cyan "  Applying $patch_name"
    if ! (cd "$target_dir" && patch -p1 < "$patch_path"); then
      print_red "Error: Failed to apply $patch_name"
      exit 1
    fi
  done
}

build_tarball() {
  local platform="$1" # mac | linux
  local out_path="$2"

  local stage_dir
  stage_dir="$(mktemp -d 2>/dev/null || mktemp -d -t basilisk-release)"

  mkdir -p "$stage_dir/basilisk"
  mkdir -p "$stage_dir/patches"

  print_cyan "Preparing $platform tarball staging dir: $stage_dir"
  copy_basilisk_source_to_stage "$stage_dir/basilisk"
  copy_selected_patches "$platform" "$stage_dir/patches"

  print_cyan "Applying patches for $platform..."
  apply_patches_in_dir "$stage_dir/basilisk" "$stage_dir/patches"

  tar -C "$stage_dir" -czf "$out_path" basilisk patches
  print_green "Created: $out_path"

  rm -rf "$stage_dir"
}

run_host_install_test() {
  local test_dir
  test_dir="$(mktemp -d 2>/dev/null || mktemp -d -t basilisk-test)"

  print_cyan "Running host install test in: $test_dir"
  cp "$REPO_ROOT/reset_install_basilisk.sh" "$test_dir/"
  cp -R "$REPO_ROOT/patches" "$test_dir/patches"
  chmod +x "$test_dir/reset_install_basilisk.sh"

  (cd "$test_dir" && ./reset_install_basilisk.sh --mode=1 --hard)
  rm -rf "$test_dir"
}

run_linux_install_test_in_docker() {
  local image_name="darcs-test"

  if ! docker info >/dev/null 2>&1; then
    print_red "Error: Docker is not running"
    exit 1
  fi

  if ! docker image inspect "$image_name" >/dev/null 2>&1; then
    print_cyan "Docker image '$image_name' not found. Building..."
    docker build -f "$REPO_ROOT/Dockerfile.darcs-test" -t "$image_name" "$REPO_ROOT"
  fi

  print_cyan "Running Linux install test in Docker..."
  docker run --rm -v "$REPO_ROOT":/repo:ro "$image_name" bash -lc "\
    set -euo pipefail && \
    mkdir -p /work && cd /work && \
    cp /repo/reset_install_basilisk.sh . && \
    cp -R /repo/patches ./patches && \
    chmod +x reset_install_basilisk.sh && \
    ./reset_install_basilisk.sh --mode=1 --hard \
  "
}

cd "$REPO_ROOT"

require_tool git
require_tool darcs
require_tool tar
require_tool patch
require_tool gh

if [[ "$SKIP_TESTS" == false ]]; then
  require_tool make
  require_tool gcc
  require_tool gawk
  require_tool curl
  if [[ "$OSTYPE" == "darwin"* ]]; then
    require_tool docker
  fi
fi

if ! gh auth status >/dev/null 2>&1; then
  print_red "Error: gh is not authenticated. Run: gh auth login"
  exit 1
fi

if [[ -n "$(git status --porcelain)" ]]; then
  print_red "Error: Working tree not clean. Commit/stash changes before releasing."
  git status --porcelain
  exit 1
fi

branch="$(git rev-parse --abbrev-ref HEAD)"
if [[ "$branch" != "main" ]]; then
  print_red "Error: Must be on 'main' branch to release (current: $branch)"
  exit 1
fi

if ! git remote get-url origin >/dev/null 2>&1; then
  print_red "Error: Missing 'origin' git remote"
  exit 1
fi

print_cyan "Fetching origin..."
git fetch origin --tags
git pull --ff-only origin main

if git ls-remote --exit-code --tags origin "refs/tags/$TAG" >/dev/null 2>&1; then
  print_red "Error: Tag already exists on origin: $TAG"
  exit 1
fi

if git rev-parse -q --verify "refs/tags/$TAG" >/dev/null; then
  if [[ "$FORCE" == true ]]; then
    print_cyan "Deleting existing local tag: $TAG"
    git tag -d "$TAG"
  else
    print_red "Error: Tag already exists locally: $TAG"
    exit 1
  fi
fi

if [[ "$DRY_RUN" == true ]]; then
  print_cyan "===== DRY RUN MODE ====="
  print_cyan "Would perform the following actions:"
  echo "  - Sync basilisk-source from upstream darcs"
  echo "  - Commit any basilisk-source changes"
  if [[ "$SKIP_TESTS" == false ]]; then
    echo "  - Run host install tests"
    if [[ "$OSTYPE" == "darwin"* ]]; then
      echo "  - Run Linux install tests in Docker"
    fi
  fi
  echo "  - Build basilisk-mac.tar.gz and basilisk-linux.tar.gz"
  echo "  - Create tag: $TAG"
  echo "  - Push to origin"
  echo "  - Create GitHub Release with assets"
  exit 0
fi

print_cyan "Syncing basilisk-source from upstream darcs..."
if [[ ! -d "$REPO_ROOT/basilisk-source/_darcs" ]]; then
  print_red "Error: Missing darcs repo at basilisk-source/_darcs"
  exit 1
fi
(cd "$REPO_ROOT/basilisk-source" && darcs pull --all)

if [[ -n "$(git status --porcelain basilisk-source)" ]]; then
  print_cyan "Committing basilisk-source changes..."
  git add basilisk-source
  git commit -m "Update basilisk-source from Darcs ($TAG)"
fi

if [[ "$SKIP_TESTS" == false ]]; then
  print_cyan "Running install tests..."
  run_host_install_test

  if [[ "$OSTYPE" == "darwin"* ]]; then
    run_linux_install_test_in_docker
  fi
fi

artifacts_dir="$(mktemp -d 2>/dev/null || mktemp -d -t basilisk-artifacts)"
print_cyan "Artifacts directory: $artifacts_dir"

mac_tar="$artifacts_dir/basilisk-mac.tar.gz"
linux_tar="$artifacts_dir/basilisk-linux.tar.gz"

build_tarball "mac" "$mac_tar"
build_tarball "linux" "$linux_tar"

sha256_file "$mac_tar" > "$mac_tar.sha256"
sha256_file "$linux_tar" > "$linux_tar.sha256"

print_green "Checksums:"
cat "$mac_tar.sha256"
cat "$linux_tar.sha256"

release_notes="$artifacts_dir/release-notes.md"
cat > "$release_notes" <<EOF
Ref-locked Basilisk release \`$TAG\`.

Assets:
- \`basilisk-mac.tar.gz\`: Basilisk source (darcs snapshot) + comphy-lab patches for macOS (local-bview not applied)
- \`basilisk-linux.tar.gz\`: Basilisk source (darcs snapshot) + comphy-lab patches for Linux (local-bview not applied)

Install (ref-locked):
- \`./reset_install_basilisk.sh --mode=4 --ref=$TAG --hard\`
- \`./reset_install_basilisk-ref-locked.sh --ref=$TAG --hard\`
EOF

print_cyan "Creating annotated tag: $TAG"
git tag -a "$TAG" -m "comphy-lab Basilisk release $TAG"

print_cyan "Pushing main and tag..."
git push origin main
git push origin "$TAG"

print_cyan "Creating GitHub Release and uploading assets..."
gh release create "$TAG" \
  "$mac_tar" "$linux_tar" "$mac_tar.sha256" "$linux_tar.sha256" \
  --title "$TAG" \
  --notes-file "$release_notes"

print_green "✅ Release published: $TAG"
print_green "Artifacts kept at: $artifacts_dir"
