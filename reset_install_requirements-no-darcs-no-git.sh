#!/bin/bash
# LEGACY WRAPPER: This script has been renamed to reset_install_basilisk-no-darcs-no-git.sh
# This wrapper exists for backwards compatibility with existing URLs.
# Please update your bookmarks/scripts to use the new name.

SCRIPT_DIR="$(cd "$(dirname "$0")" 2>/dev/null && pwd)"
NEW_SCRIPT="$SCRIPT_DIR/reset_install_basilisk-no-darcs-no-git.sh"

if [[ -f "$NEW_SCRIPT" ]]; then
    # Local execution - use the renamed script
    exec "$NEW_SCRIPT" "$@"
else
    # Remote execution (curl | bash) - fetch the new script
    exec bash <(curl -sL https://raw.githubusercontent.com/comphy-lab/basilisk-C/main/reset_install_basilisk-no-darcs-no-git.sh) "$@"
fi
