#!/usr/bin/env bash
# Shared optional-patch contract for comphy-lab Basilisk installers and releases.
#
# Default installations and release tarballs apply only ordinary platform
# patches. The two bview patches are opt-in and mutually exclusive:
#   --local-bview   -> *-local-bview.patch
#   --comphy-bview  -> *-comphy-bview.patch
#
# Source this file; do not execute it.

if [[ -n "${COMPHY_PATCH_CONTRACT_SOURCED:-}" ]]; then
  return 0 2>/dev/null || exit 0
fi
COMPHY_PATCH_CONTRACT_SOURCED=1

comphy_require_exclusive_bview_flags() {
  local apply_local="${1:-false}"
  local apply_comphy="${2:-false}"
  if [[ "$apply_local" == "true" && "$apply_comphy" == "true" ]]; then
    printf '%s\n' "Error: --local-bview and --comphy-bview are mutually exclusive." >&2
    printf '%s\n' "       --comphy-bview supersedes --local-bview; use one or the other." >&2
    return 1
  fi
  return 0
}

comphy_is_optional_bview_patch() {
  local patch_name="$1"
  [[ "$patch_name" == *"-local-bview.patch" || "$patch_name" == *"-comphy-bview.patch" ]]
}

# Print a skip reason on stdout and return 0 when the patch must not be
# applied/packaged for the requested flags. Return 1 when the patch should
# be applied.
comphy_patch_skip_reason() {
  local patch_name="$1"
  local apply_local="${2:-false}"
  local apply_comphy="${3:-false}"

  if [[ "$patch_name" == *"-local-bview.patch" && "$apply_local" != "true" ]]; then
    printf '%s\n' "use --local-bview to apply"
    return 0
  fi
  if [[ "$patch_name" == *"-comphy-bview.patch" && "$apply_comphy" != "true" ]]; then
    printf '%s\n' "use --comphy-bview to apply"
    return 0
  fi
  if [[ "$patch_name" == *"-macos-"* && "${OSTYPE:-}" != darwin* ]]; then
    printf '%s\n' "macOS-specific patch"
    return 0
  fi
  return 1
}

# Return 0 when a release tarball for platform (mac|linux) should include
# this patch file. Optional bview patches are never packaged.
comphy_release_include_patch() {
  local platform="$1"
  local patch_name="$2"

  if [[ "$patch_name" == *"-local-bview.patch" || "$patch_name" == *"-comphy-bview.patch" ]]; then
    return 1
  fi
  if [[ "$platform" == "linux" && "$patch_name" == *"-macos-"* ]]; then
    return 1
  fi
  return 0
}
