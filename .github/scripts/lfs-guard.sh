#!/usr/bin/env bash
set -euo pipefail

threshold_bytes=$((100 * 1024 * 1024))

usage() {
  cat >&2 <<'EOF'
Usage:
  .github/scripts/lfs-guard.sh fix-staged
  .github/scripts/lfs-guard.sh check-repo
EOF
}

size_bytes() {
  local file="$1"
  stat -f%z "$file" 2>/dev/null || stat -c%s "$file"
}

fix_staged() {
  if ! git lfs version >/dev/null 2>&1; then
    printf "git-lfs is not available; cannot auto-track large files.\n" >&2
    exit 1
  fi

  local tracked_new_files=0
  while IFS= read -r file; do
    [ -f "$file" ] || continue

    local size
    size=$(size_bytes "$file")
    if [ "$size" -le "$threshold_bytes" ]; then
      continue
    fi

    local attr
    attr=$(git check-attr filter -- "$file")
    if [[ "$attr" != *": filter: lfs" ]]; then
      printf "Tracking large file with Git LFS: %s\n" "$file"
      git lfs track -- "$file"
      tracked_new_files=1
    fi

    # Re-stage the file so the commit stores an LFS pointer.
    git add -- "$file"
  done < <(git diff --cached --name-only --diff-filter=AM)

  if [ "$tracked_new_files" -eq 1 ]; then
    git add .gitattributes
  fi
}

check_repo() {
  local violations=0
  while IFS=$'\t' read -r meta file; do
    set -- $meta
    local size="${4:-0}"
    case "$size" in
      ''|*[!0-9]*)
        continue
        ;;
    esac

    if [ "$size" -le "$threshold_bytes" ]; then
      continue
    fi

    local attr
    attr=$(git check-attr filter -- "$file")
    if [[ "$attr" != *": filter: lfs" ]]; then
      printf "Large file in HEAD is not tracked by Git LFS: %s (%s bytes)\n" "$file" "$size" >&2
      violations=1
    fi
  done < <(git ls-tree -r -l --full-tree HEAD)

  if [ "$violations" -eq 1 ]; then
    cat >&2 <<'EOF'

Fix:
  1) .github/scripts/setup.sh
  2) Re-add large files and commit again
EOF
    exit 1
  fi
}

mode="${1:-}"
case "$mode" in
  fix-staged)
    fix_staged
    ;;
  check-repo)
    check_repo
    ;;
  *)
    usage
    exit 2
    ;;
esac
