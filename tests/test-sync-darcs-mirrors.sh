#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
HELPER="$ROOT/scripts/sync-darcs-mirrors.sh"

fail() {
  printf 'FAIL: %s\n' "$*" >&2
  exit 1
}

pass() {
  printf 'ok %s\n' "$*"
}

run_closed_stdin() {
  "$HELPER" "$@" </dev/null
}

output="$(set +e; run_closed_stdin not-a-mirror 2>&1; echo EXIT:$?)"
printf '%s\n' "$output" | grep -q 'unknown argument' || fail "unmapped name was not rejected"
printf '%s\n' "$output" | grep -q 'EXIT:1' || fail "unmapped name did not fail"
pass "mapping rejection for unknown mirror name"

output="$(set +e; run_closed_stdin --target /tmp/evil source 2>&1; echo EXIT:$?)"
printf '%s\n' "$output" | grep -q 'arbitrary target' || fail "--target was not rejected"
printf '%s\n' "$output" | grep -q 'EXIT:1' || fail "--target did not fail"
pass "arbitrary --target is rejected"

output="$(set +e; run_closed_stdin --repo-root / source 2>&1; echo EXIT:$?)"
printf '%s\n' "$output" | grep -q 'filesystem root' || fail "repo-root / was not rejected"
printf '%s\n' "$output" | grep -q 'EXIT:1' || fail "repo-root / did not fail"
pass "filesystem-root repo-root is rejected"

missing="$(mktemp -d)"
rmdir "$missing"
output="$(set +e; run_closed_stdin --repo-root "$missing" source 2>&1; echo EXIT:$?)"
printf '%s\n' "$output" | grep -q 'not a directory' || fail "missing repo-root was not rejected"
printf '%s\n' "$output" | grep -q 'EXIT:1' || fail "missing repo-root did not fail"
pass "missing repo-root fails closed"

output="$(set +e; printf 'n\n' | "$HELPER" --target "$ROOT/../evil" wiki 2>&1; echo EXIT:$?)"
printf '%s\n' "$output" | grep -q 'arbitrary target' || fail "piped target option was accepted"
printf '%s\n' "$output" | grep -q 'EXIT:1' || fail "piped target option did not fail"
pass "no prompt or stdin dependency for rejected arguments"

symlink_root="$(mktemp -d "${TMPDIR:-/tmp}/comphy-symlink-mirror.XXXXXX")"
mkdir -p "$symlink_root/outside" "$symlink_root/repo"
ln -s "$symlink_root/outside" "$symlink_root/repo/basilisk-source"
output="$(set +e; run_closed_stdin --repo-root "$symlink_root/repo" source 2>&1; echo EXIT:$?)"
rm -rf "$symlink_root"
printf '%s\n' "$output" | grep -q 'symlink target' || fail "symlink mirror target was not rejected"
printf '%s\n' "$output" | grep -q 'EXIT:1' || fail "symlink mirror target did not fail"
pass "symlink mirror targets are rejected"

if grep -REn 'yes[[:space:]]*\|[[:space:]]*darcs' \
  "$HELPER" \
  "$ROOT/.github/workflows/sync-darcs-repositories.yml" \
  "$ROOT/release-comphy-tag.sh" >/dev/null
then
  fail "interactive yes|darcs path is still present"
fi
if grep -q 'darcs pull' "$ROOT/.github/workflows/sync-darcs-repositories.yml"; then
  fail "workflow still uses darcs pull"
fi
if grep -q 'darcs pull --all' "$ROOT/release-comphy-tag.sh"; then
  fail "release script still uses darcs pull directly"
fi
grep -q 'dont-allow-conflicts' "$HELPER" || fail "existing-repo pull is not fail-closed on conflicts"
grep -q 'size +90M' "$HELPER" || fail "helper does not reject oversized Git blobs"
grep -q 'tentative_hashed_inventory' "$HELPER" || fail "helper does not strip Darcs tentative inventory files"
pass "helper, workflow, and release script have no interactive darcs pull"

sync_one_and_check() {
  local name="$1"
  local rel_dir="$2"
  local work
  work="$(mktemp -d "$ROOT/.darcs-sync-staging.testrun.XXXXXX")"
  mkdir -p "$work/$rel_dir"
  printf 'stale\n' > "$work/$rel_dir/stale-from-test.txt"
  printf 'stale-backup\n' > "$work/$rel_dir/stale-from-test.c.~0~"

  run_closed_stdin --repo-root "$work" "$name"
  [[ -d "$work/$rel_dir/_darcs" ]] || fail "$name sync did not create a Darcs repo"
  [[ ! -e "$work/$rel_dir/stale-from-test.txt" ]] || fail "$name sync retained stale files"
  if find "$work/$rel_dir" \( -name '.*.~[0-9]*~' -o -name '*.~[0-9]*~' \) -print | grep -q .; then
    fail "$name sync introduced or retained Darcs backup files"
  fi
  [[ ! -e "$work/$rel_dir/_darcs/index" ]] || fail "$name sync left a Darcs index cache"
  [[ ! -e "$work/$rel_dir/_darcs/index.old" ]] || fail "$name sync left a Darcs index.old cache"

  tree_manifest() {
    local tree="$1"
    (
      cd "$tree" || exit 1
      find . -type f ! -path './_darcs/index' ! -path './_darcs/index.old' -print0 \
        | sort -z \
        | xargs -0 sha256sum \
        | sha256sum
    )
  }

  local first_manifest
  first_manifest="$(tree_manifest "$work/$rel_dir")"

  # Second pass must use the existing-_darcs path: discard dirt, pull, and
  # rsync the working tree without replacing the Darcs store.
  local sample
  sample="$(find "$work/$rel_dir" -type f ! -path '*/_darcs/*' -print -quit)"
  [[ -n "$sample" ]] || fail "$name sync produced no working-tree files"
  printf 'unrecorded-dirt\n' >> "$sample"
  printf 'stale-backup\n' > "$work/$rel_dir/second-pass.c.~0~"

  run_closed_stdin --repo-root "$work" "$name"
  [[ ! -e "$work/$rel_dir/second-pass.c.~0~" ]] || fail "$name incremental sync retained backup files"
  [[ ! -e "$work/$rel_dir/_darcs/patches/unrevert" ]] || fail "$name incremental sync left unrevert state"
  if find "$work/$rel_dir/_darcs" \( -name '*.tentative' -o -name 'tentative_*' \) -print | grep -q .; then
    fail "$name incremental sync left tentative Darcs state"
  fi
  if find "$work/$rel_dir" \( -name '.*.~[0-9]*~' -o -name '*.~[0-9]*~' \) -print | grep -q .; then
    fail "$name incremental sync introduced or retained Darcs backup files"
  fi
  if find "$work/$rel_dir" -type f -size +90M -print | grep -q .; then
    fail "$name incremental sync introduced a file over 90 MiB"
  fi
  local second_manifest
  second_manifest="$(tree_manifest "$work/$rel_dir")"
  [[ "$first_manifest" == "$second_manifest" ]] || fail "$name repeated sync was not idempotent"

  rm -rf "$work"
  pass "fresh $name-only sync, no backups, and repeated sync idempotence"
}

sync_existing_hashed_source() {
  local src="$ROOT/basilisk-source"
  [[ -d "$src/_darcs" ]] || fail "canonical basilisk-source is missing _darcs"
  local work
  work="$(mktemp -d "$ROOT/.darcs-sync-staging.testrun.XXXXXX")"
  mkdir -p "$work"
  rsync -a --delete \
    --exclude=_darcs/index \
    --exclude=_darcs/index.old \
    --exclude=_darcs/patches/unrevert \
    -- "$src/" "$work/basilisk-source/"
  local before_patch_count witness witness_hash
  before_patch_count="$(find "$work/basilisk-source/_darcs/patches" -type f ! -name unrevert ! -name '*.tentative' | wc -l)"
  [[ "$before_patch_count" -gt 10 ]] || fail "hashed source copy has too few patch files to be a Git-tracked store"
  witness="$(find "$work/basilisk-source/_darcs/patches" -type f ! -name unrevert ! -name '*.tentative' -print -quit)"
  [[ -n "$witness" ]] || fail "hashed source copy has no witness patch"
  witness_hash="$(sha256sum "$witness" | awk '{print $1}')"
  printf 'unrecorded-dirt\n' >> "$work/basilisk-source/src/common.h"
  printf 'stale-backup\n' > "$work/basilisk-source/hashed-pass.c.~0~"

  run_closed_stdin --repo-root "$work" source
  [[ ! -e "$work/basilisk-source/hashed-pass.c.~0~" ]] || fail "hashed incremental sync retained backup files"
  [[ ! -e "$work/basilisk-source/_darcs/patches/unrevert" ]] || fail "hashed incremental sync left unrevert state"
  if find "$work/basilisk-source/_darcs" \( -name '*.tentative' -o -name 'tentative_*' \) -print | grep -q .; then
    fail "hashed incremental sync left tentative Darcs state"
  fi
  if find "$work/basilisk-source" -type f -size +90M -print | grep -q .; then
    fail "hashed incremental sync introduced a file over 90 MiB"
  fi
  [[ -f "$witness" ]] || fail "hashed incremental sync removed a pre-existing patch file"
  [[ "$(sha256sum "$witness" | awk '{print $1}')" == "$witness_hash" ]] || fail "hashed incremental sync rewrote a pre-existing patch file"
  local after_patch_count
  after_patch_count="$(find "$work/basilisk-source/_darcs/patches" -type f ! -name unrevert ! -name '*.tentative' | wc -l)"
  [[ "$after_patch_count" -ge "$before_patch_count" ]] || fail "hashed incremental sync replaced the existing Darcs patch store"
  rm -rf "$work"
  pass "existing hashed basilisk-source stays incremental and under 90 MiB"
}

oversize="$(mktemp -d "${TMPDIR:-/tmp}/comphy-oversize-guard.XXXXXX")"
dd if=/dev/zero of="$oversize/big.bin" bs=1M count=91 status=none
find "$oversize" -type f -size +90M -print | grep -q . || fail "90 MiB size guard expression missed a 91 MiB file"
rm -rf "$oversize"
pass "size-guard find expression matches files over 90 MiB"

if [[ "${SKIP_DARCS_NETWORK_TESTS:-}" == "1" ]]; then
  pass "skipping network Darcs sync tests (SKIP_DARCS_NETWORK_TESTS=1)"
else
  sync_one_and_check source basilisk-source
  sync_existing_hashed_source
  sync_one_and_check wiki basilisk-wiki
fi
