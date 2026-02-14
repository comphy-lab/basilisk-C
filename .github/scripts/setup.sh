#!/usr/bin/env bash
set -euo pipefail

repo_root=$(git rev-parse --show-toplevel 2>/dev/null || true)
if [ -z "${repo_root}" ]; then
  printf "Run this script from inside a git working tree.\n" >&2
  exit 1
fi

cd "$repo_root"

if [ ! -d ".github/hooks" ]; then
  printf "Missing hooks directory: .github/hooks\n" >&2
  exit 1
fi

# Install Git LFS hooks/config in this user environment.
if ! git lfs version >/dev/null 2>&1; then
  printf "git-lfs is required. Install it, then re-run .github/scripts/setup.sh.\n" >&2
  exit 1
fi

git lfs install
git config --local core.hooksPath .github/hooks

configured_path=$(git config --local --get core.hooksPath)
printf "Configured local hooks path: %s\n" "$configured_path"
printf "Git LFS and hooks are now active for this clone.\n"
