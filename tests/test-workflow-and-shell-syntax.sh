#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"

fail() {
  printf 'FAIL: %s\n' "$*" >&2
  exit 1
}

pass() {
  printf 'ok %s\n' "$*"
}

workflow="$ROOT/.github/workflows/sync-darcs-repositories.yml"
[[ -f "$workflow" ]] || fail "missing workflow $workflow"

grep -q 'permissions:' "$workflow" || fail "workflow is missing permissions"
grep -q 'contents: write' "$workflow" || fail "workflow is missing contents: write"
grep -q 'concurrency:' "$workflow" || fail "workflow is missing a concurrency group"
grep -q 'group: sync-darcs-repositories' "$workflow" || fail "workflow concurrency group is wrong"
grep -q 'cancel-in-progress: false' "$workflow" || fail "workflow should not cancel an in-progress sync"
grep -q 'actions/checkout@v7' "$workflow" || fail "workflow is not using the current checkout major"
grep -q 'sync-darcs-mirrors.sh source' "$workflow" || fail "workflow does not sync the source mirror via the helper"
grep -q 'sync-darcs-mirrors.sh wiki' "$workflow" || fail "workflow does not sync the wiki mirror via the helper"
grep -q 'Auto-update basilisk-source from Darcs repository' "$workflow" || fail "workflow lost the source commit message"
grep -q 'Auto-update basilisk-wiki from Darcs repository' "$workflow" || fail "workflow lost the wiki commit message"
grep -q 'git push' "$workflow" || fail "workflow no longer pushes after separate commits"
if grep -q 'darcs pull' "$workflow"; then
  fail "workflow still uses darcs pull"
fi
pass "workflow contract (permissions, concurrency, checkout@v7, helper, commits)"

release="$ROOT/release-comphy-tag.sh"
grep -q 'scripts/comphy-patch-contract.sh' "$release" || fail "release script does not use the shared patch contract"
grep -q 'cp "$REPO_ROOT/scripts/comphy-patch-contract.sh"' "$release" || fail "host install test does not copy the patch contract"
grep -q 'cp /repo/scripts/comphy-patch-contract.sh scripts/' "$release" || fail "docker install test does not copy the patch contract"
pass "release install tests copy the shared patch contract"

if python3 -c 'import yaml' >/dev/null 2>&1; then
  python3 - "$workflow" <<'PY'
import sys
import yaml
with open(sys.argv[1], encoding="utf-8") as handle:
    data = yaml.safe_load(handle)
assert data["permissions"]["contents"] == "write"
assert data["concurrency"]["group"] == "sync-darcs-repositories"
assert data["concurrency"]["cancel-in-progress"] is False
steps = data["jobs"]["sync"]["steps"]
uses = [step.get("uses", "") for step in steps]
assert any(item.startswith("actions/checkout@v7") for item in uses), uses
PY
  pass "workflow YAML parses and matches the required contract"
else
  python3 - "$workflow" <<'PY'
from pathlib import Path
import sys
text = Path(sys.argv[1]).read_text(encoding="utf-8")
if "\t" in text:
    raise SystemExit("workflow contains tabs")
if text.count("name: Sync Darcs Repositories") != 1:
    raise SystemExit("workflow name is missing")
PY
  pass "workflow YAML basic syntax check (PyYAML not installed)"
fi

scripts=(
  "$ROOT/scripts/sync-darcs-mirrors.sh"
  "$ROOT/scripts/comphy-patch-contract.sh"
  "$ROOT/release-comphy-tag.sh"
  "$ROOT/reset_install_basilisk.sh"
  "$ROOT/reset_install_basilisk-darcs.sh"
  "$ROOT/reset_install_basilisk-no-darcs.sh"
  "$ROOT/reset_install_basilisk-no-darcs-no-git.sh"
  "$ROOT/reset_install_basilisk-ref-locked.sh"
  "$ROOT/tests/test-comphy-patch-contract.sh"
  "$ROOT/tests/test-sync-darcs-mirrors.sh"
  "$ROOT/tests/test-workflow-and-shell-syntax.sh"
  "$ROOT/tests/run.sh"
)
for script in "${scripts[@]}"; do
  bash -n "$script"
  pass "bash -n $(basename "$script")"
done
