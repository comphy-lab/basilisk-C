#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"

bash "$ROOT/test-workflow-and-shell-syntax.sh"
bash "$ROOT/test-comphy-patch-contract.sh"
bash "$ROOT/test-sync-darcs-mirrors.sh"
