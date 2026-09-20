#!/usr/bin/env bash
# Refresh the read-only basilisk.fr Darcs mirrors.
#
# A missing Git-tracked target is filled from a fresh lazy checkout. When the
# target already has `_darcs`, that hashed history is kept: replacing it from a
# lazy wiki checkout produced a pack over GitHub's 100 MiB limit. The existing
# store is aligned with upstream by discarding unrecorded dirt, obliterating
# patches that are no longer on the remote, then
# `darcs pull --all --dont-allow-conflicts`. Pull alone cannot drop extra
# patches, so a remote obliterate (wiki spam, rewritten history) previously
# left the Git mirror dirty after rsync. Unrevert state and backup files such
# as *.~0~ are not retained.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: scripts/sync-darcs-mirrors.sh [--repo-root DIR] source|wiki|all

Refresh one or both read-only upstream Darcs mirrors. Missing targets are
filled from a fresh lazy checkout. Existing Git-tracked `_darcs` histories
are updated incrementally so a lazy hashed pack cannot replace them.
Allowed mirrors:

  source  https://basilisk.fr/basilisk  ->  <repo-root>/basilisk-source
  wiki    https://basilisk.fr/wiki      ->  <repo-root>/basilisk-wiki

--repo-root defaults to this repository. Arbitrary target paths are rejected.
EOF
}

die() {
  printf 'sync-darcs-mirrors: %s\n' "$*" >&2
  exit 1
}

SCRIPT_PATH="$(cd "$(dirname "$0")" && pwd -P)/$(basename "$0")"
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_PATH")" && pwd -P)"
DEFAULT_REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd -P)"

REPO_ROOT="$DEFAULT_REPO_ROOT"
REQUESTED=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --help|-h)
      usage
      exit 0
      ;;
    --repo-root)
      [[ $# -ge 2 ]] || die "missing value for --repo-root"
      REPO_ROOT="$2"
      shift 2
      ;;
    --repo-root=*)
      REPO_ROOT="${1#*=}"
      shift
      ;;
    --target|--dest|--directory|--dir)
      die "refusing arbitrary target option '$1'; use source, wiki, or all"
      ;;
    --target=*|--dest=*|--directory=*|--dir=*)
      die "refusing arbitrary target option '${1%%=*}'"
      ;;
    source|wiki|all)
      REQUESTED+=("$1")
      shift
      ;;
    *)
      die "unknown argument: $1"
      ;;
  esac
done

if [[ ${#REQUESTED[@]} -eq 0 ]]; then
  usage >&2
  exit 1
fi

DARCS_BIN="$(command -v darcs || true)"
RSYNC_BIN="$(command -v rsync || true)"
FIND_BIN="$(command -v find || true)"
MKTEMP_BIN="$(command -v mktemp || true)"
[[ -n "$DARCS_BIN" ]] || die "darcs is required"
[[ -n "$RSYNC_BIN" ]] || die "rsync is required"
[[ -n "$FIND_BIN" ]] || die "find is required"
[[ -n "$MKTEMP_BIN" ]] || die "mktemp is required"

[[ -n "$REPO_ROOT" ]] || die "repository root is empty"
[[ -d "$REPO_ROOT" ]] || die "repository root is not a directory: $REPO_ROOT"
REPO_ROOT="$(cd "$REPO_ROOT" && pwd -P)"
[[ "$REPO_ROOT" != / ]] || die "refusing to use filesystem root as --repo-root"

MIRROR_SOURCE_URL="https://basilisk.fr/basilisk"
MIRROR_WIKI_URL="https://basilisk.fr/wiki"
MIRROR_SOURCE_DIR="basilisk-source"
MIRROR_WIKI_DIR="basilisk-wiki"

STAGING=""
cleanup_staging() {
  if [[ -n "${STAGING:-}" && -d "$STAGING" ]]; then
    rm -rf "$STAGING"
  fi
  STAGING=""
}
trap cleanup_staging EXIT INT TERM HUP

darcs_cmd() {
  # Never attach a TTY or pager. Conflicting-pull prompts must not be reachable
  # from this helper; closing stdin makes any unexpected prompt fail closed.
  env -u DARCS_EDITOR -u VISUAL -u EDITOR \
    DARCS_DONT_COLOR=1 \
    GIT_TERMINAL_PROMPT=0 \
    "$DARCS_BIN" "$@" </dev/null
}

resolve_mirror() {
  local name="$1"
  case "$name" in
    source)
      printf '%s\t%s\n' "$MIRROR_SOURCE_URL" "$MIRROR_SOURCE_DIR"
      ;;
    wiki)
      printf '%s\t%s\n' "$MIRROR_WIKI_URL" "$MIRROR_WIKI_DIR"
      ;;
    *)
      die "refusing unmapped mirror name: $name"
      ;;
  esac
}

assert_whitelisted_target() {
  local target="$1"
  local expected_source="$REPO_ROOT/$MIRROR_SOURCE_DIR"
  local expected_wiki="$REPO_ROOT/$MIRROR_WIKI_DIR"
  local physical
  if [[ "$target" != "$expected_source" && "$target" != "$expected_wiki" ]]; then
    die "refusing non-whitelisted target: $target"
  fi
  if [[ -L "$target" ]]; then
    die "refusing symlink target: $target"
  fi
  if [[ -e "$target" ]]; then
    physical="$(cd "$target" && pwd -P)"
    if [[ "$physical" != "$target" ]]; then
      die "refusing target that resolves outside the whitelist: $physical"
    fi
  fi
}

darcs_is_clean() {
  local repo="$1"
  local output
  output="$(darcs_cmd whatsnew --repodir "$repo" --summary --look-for-adds 2>&1)" || true
  if [[ -z "${output//[$'\t\r\n ']}" ]]; then
    return 0
  fi
  if printf '%s\n' "$output" | grep -Eq '^No changes!?$'; then
    return 0
  fi
  printf 'sync-darcs-mirrors: checkout is not Darcs-clean:\n%s\n' "$output" >&2
  return 1
}

assert_no_darcs_backups() {
  local root="$1"
  local found
  found="$("$FIND_BIN" "$root" \( -name '.*.~[0-9]*~' -o -name '*.~[0-9]*~' \) -print)"
  if [[ -n "$found" ]]; then
    printf 'sync-darcs-mirrors: refusing Darcs backup debris:\n%s\n' "$found" >&2
    return 1
  fi
  return 0
}

git_lfs_attr_matches() {
  local relpath="$1"
  local attrfile="$REPO_ROOT/.gitattributes"
  local output value line pattern rest
  if command -v git >/dev/null 2>&1 \
    && git -C "$REPO_ROOT" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    output="$(git -C "$REPO_ROOT" check-attr filter -- "$relpath" 2>/dev/null || true)"
    value="${output##*: }"
    [[ "$value" == lfs ]] && return 0
    return 1
  fi
  [[ -f "$attrfile" ]] || return 1
  value=""
  while IFS= read -r line || [[ -n "$line" ]]; do
    line="${line%%#*}"
    line="${line#"${line%%[![:space:]]*}"}"
    [[ -n "$line" ]] || continue
    pattern="${line%%[[:space:]]*}"
    rest="${line#"$pattern"}"
    case "$relpath" in
      $pattern)
        case "$rest" in
          *-filter*) value="" ;;
          *filter=lfs*) value=lfs ;;
          *filter=*) value="" ;;
        esac
        ;;
    esac
  done < "$attrfile"
  [[ "$value" == lfs ]]
}

assert_no_oversized_git_files() {
  local root="$1"
  local found="" path relpath
  # GitHub rejects ordinary blobs over 100 MiB. Files already declared as
  # Git LFS in .gitattributes are stored as pointers and are allowed.
  while IFS= read -r -d '' path; do
    relpath="${path#"$REPO_ROOT"/}"
    if git_lfs_attr_matches "$relpath"; then
      printf 'sync-darcs-mirrors: allowing Git LFS file over 90 MiB: %s\n' "$relpath"
      continue
    fi
    found+="$path"$'\n'
  done < <("$FIND_BIN" "$root" -type f -size +90M -print0)
  if [[ -n "$found" ]]; then
    printf 'sync-darcs-mirrors: refusing files over 90 MiB:\n%s' "$found" >&2
    return 1
  fi
  return 0
}

normalize_darcs_cache() {
  local repo="$1"
  # Rebuildable working-tree caches and revert/rebase undo state. Leaving
  # them in Git produces host-specific churn after every sync.
  rm -f \
    "$repo/_darcs/index" \
    "$repo/_darcs/index.old" \
    "$repo/_darcs/lock" \
    "$repo/_darcs/rebase.tentative" \
    "$repo/_darcs/tentative_hashed_inventory" \
    "$repo/_darcs/tentative_pristine"
  discard_darcs_unrevert "$repo"
  if [[ -d "$repo/_darcs/patches" ]]; then
    "$FIND_BIN" "$repo/_darcs/patches" -maxdepth 1 -type f \
      -name '*.tentative' -delete
  fi
}

discard_darcs_unrevert() {
  local repo="$1"
  # Revert writes an unrevert bundle. Obliterate then prompts
  # "This operation will make unrevert impossible!" even with --all.
  rm -f "$repo/_darcs/patches/unrevert"
}

remove_darcs_backups() {
  local root="$1"
  [[ -d "$root" ]] || return 0
  # Revert backups can be files or directories (name.~0~). Delete both.
  "$FIND_BIN" "$root" -depth \
    \( -name '.*.~[0-9]*~' -o -name '*.~[0-9]*~' \) \
    -exec rm -rf -- {} +
}

finish_mirror() {
  local name="$1"
  local url="$2"
  local rel_dir="$3"
  local target="$4"

  [[ -d "$target/_darcs" ]] || die "replaced $name mirror is missing _darcs"
  darcs_is_clean "$target" || die "replaced $name mirror is not Darcs-clean"
  assert_no_darcs_backups "$target" || die "replaced $name mirror contained backup files"
  # `darcs whatsnew` recreates rebuildable caches during validation.
  # Strip them before the size scan so they cannot enter Git.
  normalize_darcs_cache "$target"
  assert_no_oversized_git_files "$target" || die "replaced $name mirror contained oversized files"
  cleanup_staging
  printf 'sync-darcs-mirrors: %s is a clean read-only mirror of %s\n' "$rel_dir" "$url"
}

sync_existing_hashed() {
  local name="$1"
  local url="$2"
  local target="$3"

  [[ ! -L "$target/_darcs" ]] || die "refusing symlink $target/_darcs"
  assert_whitelisted_target "$target"
  printf 'sync-darcs-mirrors: discarding unrecorded changes in %s\n' "$target"
  # A clean tree exits 0 with nothing to revert; a dirty tree must not prompt.
  if ! darcs_cmd revert --repodir "$target" --all; then
    die "darcs revert failed for $target"
  fi
  remove_darcs_backups "$target"
  discard_darcs_unrevert "$target"
  printf 'sync-darcs-mirrors: dropping patches not in %s from %s\n' "$url" "$target"
  if ! darcs_cmd obliterate --repodir "$target" --not-in-remote="$url" --all; then
    die "darcs obliterate of patches not on $url failed for $target"
  fi
  printf 'sync-darcs-mirrors: pulling %s into existing %s\n' "$url" "$target"
  if ! darcs_cmd pull --repodir "$target" --all --dont-allow-conflicts "$url"; then
    die "darcs pull failed for $name ($url)"
  fi
  if ! darcs_cmd revert --repodir "$target" --all; then
    die "darcs revert after pull failed for $target"
  fi
  remove_darcs_backups "$target"
  discard_darcs_unrevert "$target"
  printf 'sync-darcs-mirrors: garbage-collecting unreferenced hashed files in %s\n' "$target"
  if ! darcs_cmd optimize clean --repodir "$target"; then
    die "darcs optimize clean failed for $target"
  fi
}

sync_fresh_lazy() {
  local name="$1"
  local url="$2"
  local target="$3"
  local checkout

  STAGING="$("$MKTEMP_BIN" -d "$REPO_ROOT/.darcs-sync-staging.XXXXXX")"
  checkout="$STAGING/checkout"

  printf 'sync-darcs-mirrors: lazy checkout %s -> %s\n' "$url" "$checkout"
  if ! darcs_cmd get --lazy "$url" "$checkout"; then
    die "darcs get failed for $name ($url)"
  fi
  [[ -d "$checkout/_darcs" ]] || die "fresh $name checkout is missing _darcs"
  darcs_is_clean "$checkout" || die "fresh $name checkout was not Darcs-clean"
  assert_no_darcs_backups "$checkout" || die "fresh $name checkout contained backup files"

  mkdir -p "$target"
  assert_whitelisted_target "$target"
  printf 'sync-darcs-mirrors: replacing %s\n' "$target"
  if ! "$RSYNC_BIN" -a --delete -- "$checkout"/ "$target"/; then
    die "rsync replace failed for $target"
  fi
}

sync_one() {
  local name="$1"
  local mapping url rel_dir target

  mapping="$(resolve_mirror "$name")"
  url="${mapping%%$'\t'*}"
  rel_dir="${mapping#*$'\t'}"
  target="$REPO_ROOT/$rel_dir"
  assert_whitelisted_target "$target"

  if [[ -d "$target/_darcs" ]]; then
    sync_existing_hashed "$name" "$url" "$target"
  else
    sync_fresh_lazy "$name" "$url" "$target"
  fi

  finish_mirror "$name" "$url" "$rel_dir" "$target"
}

MIRRORS=()
for item in "${REQUESTED[@]}"; do
  case "$item" in
    all)
      MIRRORS+=(source wiki)
      ;;
    source|wiki)
      MIRRORS+=("$item")
      ;;
    *)
      die "refusing unmapped mirror name: $item"
      ;;
  esac
done

# Deduplicate while preserving order.
SEEN_SOURCE=0
SEEN_WIKI=0
for item in "${MIRRORS[@]}"; do
  case "$item" in
    source)
      [[ "$SEEN_SOURCE" -eq 0 ]] || continue
      SEEN_SOURCE=1
      sync_one source
      ;;
    wiki)
      [[ "$SEEN_WIKI" -eq 0 ]] || continue
      SEEN_WIKI=1
      sync_one wiki
      ;;
  esac
done
