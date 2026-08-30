#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck source=../scripts/comphy-patch-contract.sh
. "$ROOT/scripts/comphy-patch-contract.sh"

fail() {
  printf 'FAIL: %s\n' "$*" >&2
  exit 1
}

pass() {
  printf 'ok %s\n' "$*"
}

selected_patches() {
  local apply_local="${1:-false}"
  local apply_comphy="${2:-false}"
  local patch_file patch_name
  for patch_file in "$ROOT"/patches/*.patch; do
    patch_name="$(basename "$patch_file")"
    if comphy_patch_skip_reason "$patch_name" "$apply_local" "$apply_comphy" >/dev/null; then
      continue
    fi
    printf '%s\n' "$patch_name"
  done
}

release_patches() {
  local platform="$1"
  local patch_file patch_name
  for patch_file in "$ROOT"/patches/*.patch; do
    patch_name="$(basename "$patch_file")"
    if comphy_release_include_patch "$platform" "$patch_name"; then
      printf '%s\n' "$patch_name"
    fi
  done
}

assert_contains() {
  local haystack="$1"
  local needle="$2"
  printf '%s\n' "$haystack" | grep -Fxq "$needle" || fail "expected '$needle' in selection"
}

assert_missing() {
  local haystack="$1"
  local needle="$2"
  if printf '%s\n' "$haystack" | grep -Fxq "$needle"; then
    fail "did not expect '$needle' in selection"
  fi
}

default_sel="$(selected_patches false false)"
assert_missing "$default_sel" "2026-01-06-local-bview.patch"
assert_missing "$default_sel" "2026-08-14-comphy-bview.patch"
if [[ "${OSTYPE:-}" == darwin* ]]; then
  assert_contains "$default_sel" "2025-11-03-macos-mman-compatibility.patch"
fi
pass "default patch selection excludes both optional bview patches"

local_sel="$(selected_patches true false)"
assert_contains "$local_sel" "2026-01-06-local-bview.patch"
assert_missing "$local_sel" "2026-08-14-comphy-bview.patch"
pass "default plus --local-bview selects only local-bview"

comphy_sel="$(selected_patches false true)"
assert_contains "$comphy_sel" "2026-08-14-comphy-bview.patch"
assert_missing "$comphy_sel" "2026-01-06-local-bview.patch"
pass "default plus --comphy-bview selects only comphy-bview"

if comphy_require_exclusive_bview_flags true true >/dev/null 2>&1; then
  fail "both optional flags should fail"
fi
pass "passing both optional flags fails early"

for installer in \
  reset_install_basilisk.sh \
  reset_install_basilisk-darcs.sh \
  reset_install_basilisk-no-darcs.sh \
  reset_install_basilisk-no-darcs-no-git.sh \
  reset_install_basilisk-ref-locked.sh
do
  output="$(set +e; "$ROOT/$installer" --local-bview --comphy-bview 2>&1; echo EXIT:$?)"
  printf '%s\n' "$output" | grep -q 'mutually exclusive' || fail "$installer did not report mutual exclusion"
  printf '%s\n' "$output" | grep -q 'EXIT:1' || fail "$installer did not fail when both optional flags were set"
  pass "$installer rejects both optional flags"
done

standalone="$(mktemp -d "${TMPDIR:-/tmp}/comphy-standalone-installer.XXXXXX")"
packaging="$(mktemp -d "${TMPDIR:-/tmp}/comphy-release-packaging.XXXXXX")"
trap 'rm -rf "$packaging" "$standalone"' EXIT
for installer in \
  reset_install_basilisk.sh \
  reset_install_basilisk-darcs.sh \
  reset_install_basilisk-no-darcs.sh \
  reset_install_basilisk-no-darcs-no-git.sh
do
  cp "$ROOT/$installer" "$standalone/"
  output="$(set +e; "$standalone/$installer" --hard 2>&1; echo EXIT:$?)"
  printf '%s\n' "$output" | grep -q 'comphy-patch-contract.sh' || fail "$installer did not require the shared contract"
  printf '%s\n' "$output" | grep -q 'EXIT:1' || fail "$installer did not fail closed without the shared contract"
  pass "$installer fails closed without the shared contract"
done
help_output="$(set +e; "$standalone/reset_install_basilisk.sh" --help 2>&1; echo EXIT:$?)"
printf '%s\n' "$help_output" | grep -q 'Usage:' || fail "main installer --help failed without the shared contract"
printf '%s\n' "$help_output" | grep -q 'EXIT:0' || fail "main installer --help did not succeed without the shared contract"
pass "main installer --help works without the shared contract"

for platform in mac linux; do
  packaged="$(release_patches "$platform")"
  assert_missing "$packaged" "2026-01-06-local-bview.patch"
  assert_missing "$packaged" "2026-08-14-comphy-bview.patch"
done
assert_contains "$(release_patches mac)" "2025-11-03-macos-mman-compatibility.patch"
assert_missing "$(release_patches linux)" "2025-11-03-macos-mman-compatibility.patch"
pass "default release tarball selection excludes both optional bview patches"

for platform in mac linux; do
  dest="$packaging/$platform"
  mkdir -p "$dest"
  for patch_file in "$ROOT"/patches/*.patch; do
    patch_name="$(basename "$patch_file")"
    if comphy_release_include_patch "$platform" "$patch_name"; then
      cp "$patch_file" "$dest/"
    fi
  done
  [[ ! -e "$dest/2026-01-06-local-bview.patch" ]] || fail "$platform tarball staged local-bview"
  [[ ! -e "$dest/2026-08-14-comphy-bview.patch" ]] || fail "$platform tarball staged comphy-bview"
  if [[ "$platform" == linux ]]; then
    [[ ! -e "$dest/2025-11-03-macos-mman-compatibility.patch" ]] || fail "linux tarball staged macos patch"
    [[ -z "$(ls -A "$dest")" ]] || fail "linux default tarball should stage no patches"
  else
    [[ -f "$dest/2025-11-03-macos-mman-compatibility.patch" ]] || fail "mac tarball missing macos patch"
  fi
  pass "$platform release staging copies only the default platform patches"
done
