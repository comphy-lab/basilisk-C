# Basilisk Patches

This directory contains patches maintained by comphy-lab that fix bugs or add features not yet available in upstream Basilisk. These patches are automatically applied by the installation scripts (`reset_install_basilisk*.sh`).

## Patch Naming Convention

Patches follow the format: `YYYY-MM-DD-descriptive-name.patch`

The date indicates when the patch was created. Patches are applied in alphabetical order (which corresponds to chronological order due to the naming convention).

## Patch Format

Each patch includes a metadata header followed by standard unified diff format:

```diff
# Patch: feature-description
# Author: Name <email@example.com>
# Date: YYYY-MM-DD
# Description: Brief description of what this patch does
#
# Optional longer explanation...
#
# Platform: All | macOS | Linux
# Files: list of modified files
#
diff -rN -u old-basilisk/path/to/file new-basilisk/path/to/file
--- old-basilisk/path/to/file    timestamp
+++ new-basilisk/path/to/file    timestamp
@@ -line,count +line,count @@
...
```

Lines starting with `#` are comments, ignored by `patch` and `darcs apply`. This format provides traceability while remaining compatible with standard patching tools.

## Available Patches

### 2025-11-03-macos-mman-compatibility.patch

**Platform:** macOS only
**Files modified:** `src/ast/std/sys/mman.h`
**Credit:** Peter Croxford (issue identification)

Fixes compilation errors on macOS related to memory mapping. The macOS SDK doesn't expose certain POSIX memory mapping constants and function declarations that Basilisk expects.

**What it fixes:**
- Defines `MAP_ANON` and `MAP_ANONYMOUS` constants for anonymous memory mapping
- Defines `POSIX_MADV_DONTNEED` and `MADV_DONTNEED` for memory advice
- Explicitly declares the `madvise()` function

**Symptoms without patch:**
```
error: use of undeclared identifier 'MAP_ANONYMOUS'
error: implicit declaration of function 'madvise'
```

---

### 2026-01-06-local-bview.patch

**Platform:** All
**Files modified:** `src/bview.c`, `src/display.h`
**Related:** [bview-local-client](https://github.com/comphy-lab/bview-local-client)

Adds support for using a locally-hosted JavaScript client with bview instead of the remote client at basilisk.fr. This enables offline visualization and faster loading times during development.

**What it adds:**
- `--local` flag: Uses `localhost:8000` for the JavaScript client
- `--local=PORT` flag: Uses `localhost:PORT` for custom port configurations
- `display_set_js_base()` function for runtime URL override
- Proper usage help with `-h` flag

**Usage:**
```shell
# Start local client server first
cd bview-local-client && python3 -m http.server 8000

# Then run bview with --local
bview2D --local dump
bview3D --local=3000 dump
```

**Note:** The `--local` flag only changes where the JavaScript client is served from. The WebSocket connection always connects back to your local bview process.

---

### 2026-01-13-mpi-tree-dump-header-fix.patch

**Platform:** All (affects MPI builds)
**Files modified:** `src/output.h`

Fixes uninitialized `header.n` values in the tree dump functionality. Previously, the grid dimensions were only set inside the `#if MULTIGRID_MPI` block, leaving them uninitialized for non-MPI builds or causing incorrect behavior.

**What it fixes:**
```c
// Before (incorrect): header.n only set inside MULTIGRID_MPI block
#if MULTIGRID_MPI
  foreach_dimension()
    header.n.x = Dimensions.x;
  MPI_Barrier (MPI_COMM_WORLD);
#endif

// After (correct): header.n set unconditionally
  foreach_dimension()
    header.n.x = Dimensions.x;
#if MULTIGRID_MPI
  MPI_Barrier (MPI_COMM_WORLD);
#endif
```

**Symptoms without patch:**
- Incorrect or garbage values in dump file headers
- Potential issues when restoring dump files

---

## Applying Patches Manually

To apply a patch to an existing Basilisk installation:

```shell
cd /path/to/basilisk
curl -O https://raw.githubusercontent.com/comphy-lab/basilisk-C/main/patches/PATCH_NAME.patch
patch -p1 < PATCH_NAME.patch

# Rebuild
cd src
make clean
make
```

To reverse a patch:
```shell
patch -R -p1 < PATCH_NAME.patch
```

## Creating New Patches

When creating a new patch for this repository:

1. Set up temporary directories:
   ```shell
   PATCH_WORK=$(mktemp -d)
   mkdir -p "$PATCH_WORK/old-basilisk/src" "$PATCH_WORK/new-basilisk/src"
   ```

2. Copy source files to both directories:
   ```shell
   cp /path/to/basilisk/src/file.c "$PATCH_WORK/old-basilisk/src/"
   cp /path/to/basilisk/src/file.c "$PATCH_WORK/new-basilisk/src/"
   ```

3. Make changes only to files in `new-basilisk/`

4. Create metadata header and generate patch:
   ```shell
   cat > "$PATCH_WORK/header.txt" << 'EOF'
   # Patch: your-patch-name
   # Author: Your Name <your@email.com>
   # Date: YYYY-MM-DD
   # Description: Brief description
   #
   # Platform: All
   # Files: src/file.c
   #
   EOF

   cd "$PATCH_WORK"
   cat header.txt > patch.diff
   diff -rN -u old-basilisk new-basilisk >> patch.diff
   ```

5. Test the patch applies cleanly:
   ```shell
   patch --dry-run -p1 -d /path/to/basilisk < patch.diff
   ```

6. Save to patches directory and update this README:
   ```shell
   cp patch.diff /path/to/patches/YYYY-MM-DD-description.patch
   ```

## Upstream Status

| Patch | Submitted Upstream | Status |
|-------|-------------------|--------|
| macos-mman-compatibility | No | comphy-lab maintained |
| local-bview | No | comphy-lab maintained |
| mpi-tree-dump-header-fix | No | comphy-lab maintained |

Patches marked as "comphy-lab maintained" are specific to our workflow or haven't been submitted upstream yet. If you believe a patch should be submitted to the main Basilisk project, please open an issue.
