/**
# `multigrid.h` / `multigrid3D.h` MPI crashes in `boundary_stencil()` when `max(nx, ny, nz) >= 3` and `uf.n[wall] = 0`

On a 2D or 3D MPI multigrid (`grid/multigrid.h` or `grid/multigrid3D.h`),
imposing normal face-velocity boundary conditions (such as
`uf.n[left] = 0; uf.n[right] = 0;`) crashes immediately with a `PMPI_Recv`
error or glibc `sysmalloc` heap corruption whenever any axis is split
across **3 or more MPI ranks** (`max(nx, ny, nz) >= 3`).

With at most 2 ranks per axis (`1x2x1`, `2x2x1`, `2x2x2`, etc.), or when
the `uf.n` boundary conditions are omitted (`-DNO_UF=1`), the same code
runs cleanly.

## How to reproduce

~~~bash
# Fails whenever max(nx, ny, nz) >= 3:
CC='mpicc -D_MPI=3' CFLAGS='-O2 -DNX=1 -DNY=3 -DNZ=1' make multigrid_mpi_uf_bc.tst
CC='mpicc -D_MPI=4' CFLAGS='-O2 -DNX=1 -DNY=4 -DNZ=1' make multigrid_mpi_uf_bc.tst
CC='mpicc -D_MPI=4' CFLAGS='-O2 -DUSE_2D=1 -DNX=4 -DNY=1' make multigrid_mpi_uf_bc.tst
CC='mpicc -D_MPI=8' CFLAGS='-O2 -DNX=2 -DNY=4 -DNZ=1' make multigrid_mpi_uf_bc.tst

# Passes when max(nx, ny, nz) <= 2 or when uf.n[wall] BCs are omitted:
CC='mpicc -D_MPI=8' CFLAGS='-O2 -DNX=2 -DNY=2 -DNZ=2' make multigrid_mpi_uf_bc.tst
CC='mpicc -D_MPI=4' CFLAGS='-O2 -DNX=1 -DNY=4 -DNZ=1 -DNO_UF=1' make multigrid_mpi_uf_bc.tst
~~~

Verified configurations:

| Config | `np` | `max(nx, ny, nz)` | `uf.n[wall] = 0` | Result |
| :--- | :-: | :-: | :-: | :---: |
| `1x2x1`, `2x1x1`, `1x1x2`, `2x2x1`, `2x2x2` | 2, 4, 8 | 2 | Yes | PASS |
| `1x3x1` | 3 | 3 | Yes | FAIL |
| `1x4x1`, `4x1x1`, `1x1x4` | 4 | 4 | Yes | FAIL |
| `2x4x1`, `4x2x1`, `1x2x4` | 8 | 4 | Yes | FAIL |
| 2D `4x1` (`grid/multigrid.h`) | 4 | 4 | Yes | FAIL |
| 3D `1x4x1` & 2D `4x1` (`-DNO_UF=1`) | 4 | 4 | No | PASS |
*/

#if USE_2D
# include "grid/multigrid.h"
#else
# include "grid/multigrid3D.h"
#endif
#include "navier-stokes/centered.h"

#ifndef NX
# define NX 1
#endif
#ifndef NY
# define NY 4
#endif
#ifndef NZ
# define NZ 1
#endif

#if !NO_UF
uf.n[left]   = 0; uf.n[right]  = 0;
uf.n[bottom] = 0; uf.n[top]    = 0;
# if dimension == 3
uf.n[back]   = 0; uf.n[front]  = 0;
# endif
#endif

int main()
{
#if dimension == 2
  dimensions (nx = NX, ny = NY);
#else
  dimensions (nx = NX, ny = NY, nz = NZ);
#endif
  init_grid (NX*(1 << 4));
  run();
}

event end (i = 2)
{
  if (pid() == 0)
    fprintf (stderr, "PASS [%dx%dx%d, np=%d]: N=%d tn=%ld\n",
             NX, NY, NZ, npe(), N, grid->tn);
}

/**
## Root cause

In [/src/grid/stencils.h](/src/grid/stencils.h) (lines 122–138),
`boundary_stencil()` decides for each dirty face-vector component whether
to exchange it via **full** `boundary()` (`loop->listc`) or via
**face-only** `boundary_face()` (`loop->listf`). It makes that decision
locally on each MPI rank by testing whether either opposite face along
that axis has a non-NULL normal boundary condition:

~~~c
scalar sn = s.v.x.i >= 0 ? s.v.x : s;
foreach_dimension()
  if (s.v.x.i == s.i) {
    if (sn.boundary[left] || sn.boundary[right])
      loop->listc = list_append (loop->listc, s);
    else if (!(s.stencil.bc & s_face))
      loop->listf.x = list_append (loop->listf.x, s);
  }
~~~

(In Basilisk, the normal-component BC of a vector `v` on any face `d` is
stored in `v.x.boundary[d]`.)

When the user sets `uf.n[left] = 0; uf.n[right] = 0;` (etc.), every rank
starts with non-NULL normal BCs on all walls. During `init_grid()`,
`mpi_boundary_new()` in [/src/grid/multigrid-mpi.h](/src/grid/multigrid-mpi.h)
(lines 201–208) calls `periodic_boundary(d)` in
[/src/grid/cartesian-common.h](/src/grid/cartesian-common.h) (lines
1049–1053) on every **internal MPI interface**, setting
`v.x.boundary[d] = NULL` there so physical wall BCs are not applied
between neighboring ranks.

Now consider 4 ranks along an axis (e.g. `dimensions(4, 1, 1)`):

* **Wall ranks (`0` and `3`)** still touch one physical wall, so one of
  `sn.boundary[left]` or `sn.boundary[right]` stays non-NULL. They append
  `uf.x` to `loop->listc` and exchange it in `boundary(listc)`.
* **Interior ranks (`1` and `2`)** have MPI neighbors on both sides, so
  *both* `sn.boundary[left]` and `sn.boundary[right]` are cleared to
  `NULL`. They append `uf.x` to `loop->listf` and exchange it in
  `boundary_face(listf)` first.

Because neighboring ranks place `uf` in different lists, they enter MPI
halo exchanges with mismatched scalar counts (`snd_*` / `rcv_*` in
[/src/grid/multigrid-mpi.h](/src/grid/multigrid-mpi.h), line 114),
causing an immediate `PMPI_Recv` error or glibc `sysmalloc` crash.
With at most 2 ranks per axis, every rank touches a physical wall, so all
ranks agree and the bug stays hidden.

## Proposed fix

Synchronize the `sn.boundary[left] || sn.boundary[right]` decision across
MPI ranks on non-tree multigrids in
[/src/grid/stencils.h](/src/grid/stencils.h) so that if any rank requires
full BCs on a face component, all ranks append it to `loop->listc`:

~~~diff
--- a/src/grid/stencils.h
+++ b/src/grid/stencils.h
@@ -133,4 +133,8 @@
-		if (sn.boundary[left] || sn.boundary[right])
+		int has_bc = (sn.boundary[left] || sn.boundary[right]);
+#if _MPI && !TREE
+		mpi_all_reduce (has_bc, MPI_INT, MPI_MAX);
+#endif
+		if (has_bc)
 		  loop->listc = list_append (loop->listc, s);
 		else if (!(s.stencil.bc & s_face))
 		  loop->listf.x = list_append (loop->listf.x, s);
~~~

(Alternatively, as noted by the existing `fixme` above line 133 of
`stencils.h`, imposing flux BCs directly inside `boundary_face()` would
avoid splitting face vectors between `listc` and `listf` altogether.)
*/
