#!/bin/bash
set -euo pipefail

hal_id="${1:-}"

case "$hal_id" in
    ""|*[!A-Za-z0-9._:-]*)
        printf '@misc{error, title={invalid HAL id}}\n'
        exit 0
        ;;
esac

# cache="$HOME/.halcache"
# if ! found=$(grep -F -- "$hal_id" "$cache" 2> /dev/null); then
    printf 'looking for %s on HAL...\n' "$hal_id" >&2
    url="https://api.archives-ouvertes.fr/search/?q=halId_s:${hal_id}&wt=bibtex&rows=1"
    found=$(wget --user-agent="" -q -O - "$url" || true)
    if [[ -z "$found" ]]; then
        found="@misc{error, title={${hal_id} not found or server error}}"
        printf '%s\n' "$found" >&2
#    else
#	echo "$found" >> "$cache"
    fi
# fi
printf '%s\n' "$found"
