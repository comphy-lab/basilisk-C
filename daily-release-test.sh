#!/bin/bash
# Daily automated release test runner
# Runs release-comphy-tag.sh --test-only and opens a GitHub issue on failure

set -uo pipefail

# Ensure Homebrew paths are available (cron has minimal PATH)
if [[ -d "/opt/homebrew/bin" ]]; then
  export PATH="/opt/homebrew/bin:/opt/homebrew/sbin:$PATH"
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
LOGS_DIR="$SCRIPT_DIR/logs"
DATE_STR="$(date -u +%Y-%m-%d)"
TIME_STR="$(date -u +%H:%M:%S)"
LOG_FILE="$LOGS_DIR/release-test-$DATE_STR.log"

mkdir -p "$LOGS_DIR"

log() {
  echo "[$(date -u +%Y-%m-%d\ %H:%M:%S)] $1" | tee -a "$LOG_FILE"
}

cleanup_old_logs() {
  # Keep only last 30 days of logs
  find "$LOGS_DIR" -name "release-test-*.log" -mtime +30 -delete 2>/dev/null || true
}

open_issue_on_failure() {
  local exit_code="$1"

  # Check if an open issue with this label already exists
  local existing_issue
  existing_issue=$(gh issue list --repo comphy-lab/basilisk-C --label "automated-test-failure" --state open --limit 1 --json number --jq '.[0].number // empty' 2>/dev/null || true)

  if [[ -n "$existing_issue" ]]; then
    log "Open issue #$existing_issue already exists, adding comment instead of creating new issue"
    local comment_body="**Test failed again on $DATE_STR at $TIME_STR UTC**

Exit code: \`$exit_code\`

<details>
<summary>Last 80 lines of log</summary>

\`\`\`
$(tail -80 "$LOG_FILE")
\`\`\`

</details>"
    gh issue comment "$existing_issue" --repo comphy-lab/basilisk-C --body "$comment_body"
    return
  fi

  # Create new issue
  local issue_body="The daily automated release test failed on **$DATE_STR** at **$TIME_STR UTC**.

**Exit code:** \`$exit_code\`

## Log output (last 80 lines)

\`\`\`
$(tail -80 "$LOG_FILE")
\`\`\`

---
*This issue was created automatically by \`daily-release-test.sh\`*"

  log "Creating GitHub issue for test failure..."
  gh issue create \
    --repo comphy-lab/basilisk-C \
    --title "[Automated] Release test failed on $DATE_STR" \
    --body "$issue_body" \
    --label "automated-test-failure"
}

# Main execution
log "=========================================="
log "Starting daily release test"
log "=========================================="

cleanup_old_logs

cd "$SCRIPT_DIR"

# Run the test and capture output (use tee for verbose console output)
log "Running: ./release-comphy-tag.sh --test-only"
if ./release-comphy-tag.sh --test-only 2>&1 | tee -a "$LOG_FILE"; then
  log "SUCCESS: Release test passed"
  exit 0
else
  exit_code=$?
  log "FAILURE: Release test failed with exit code $exit_code"
  open_issue_on_failure "$exit_code"
  exit $exit_code
fi
