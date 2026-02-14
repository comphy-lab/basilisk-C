---
name: basilisk-release-tag
description: Cut and publish a new Basilisk release tag in this repository using `release-comphy-tag.sh`. Use when asked to make a release, cut a tag, publish release assets, or run release tests for this repo.
---

# Basilisk Release Tag

Run this workflow from the repository root.

## Commands

- Default release (UTC date tag): `./release-comphy-tag.sh`
- Explicit tag: `./release-comphy-tag.sh --tag=vYYYY-MM-DD`
- Preflight only: `./release-comphy-tag.sh --dry-run`
- Test only (no tag/release): `./release-comphy-tag.sh --test-only`

## Workflow

1. Verify release preconditions:
   - `git branch --show-current` is `main`
   - `git status --porcelain` is empty
   - `gh auth status` succeeds
2. Optionally run a preflight:
   - `./release-comphy-tag.sh --dry-run`
3. Run the release script with default tag or `--tag=...`.
4. Confirm outcome and report:
   - published tag
   - release URL
   - uploaded assets: `basilisk-mac.tar.gz`, `basilisk-linux.tar.gz`, and matching `.sha256` files
5. Verify release details:
   - `gh release view <tag> --repo comphy-lab/basilisk-C`

## Guardrails

- Do not pass `--skip-tests` unless explicitly requested.
- Do not pass `--force` unless explicitly requested.
- Use `release-comphy-tag.sh` as the single release path; do not replace with ad-hoc git/gh commands.
