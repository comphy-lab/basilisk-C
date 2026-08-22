/**
# restore() fails on multigrid when dimensions() is not a power of two

When the domain aspect ratio set with *dimensions()* is not a power of
two, *dump()* followed by *restore()* fails with

```
grid depths do not match. Aborting.
```

even though the dump and the restore use the same executable and the
same parameters. Power-of-two layouts (1:1, 2:1, 4:1, ...) are not
affected, which is probably why this went unnoticed.

*/

#include "grid/multigrid.h"
#include "run.h"

scalar s[];

int main()
{
  dimensions (nx = 3, ny = 2); // 2:1 or 4:1 instead and the bug disappears
  init_grid (64);
  fprintf (stderr, "after init_grid: depth = %d, N = %d\n", depth(), N);

  foreach()
    s[] = x + y;

  dump ("dump");
  restore (file = "dump");
  fprintf (stderr, "restored:        depth = %d, N = %d\n", depth(), N);
}

/**

The issue can be traced to output.h. *init_grid()* derives the depth
from its argument with a rounding-up integer log,

~~~c
grid->depth = log_base2 (n/Dimensions.x);
N = (1 << grid->depth)*Dimensions.x;
~~~

so `dimensions (nx = 3, ny = 2)` with `init_grid (64)` gives depth 5,
i.e. a 96 x 64 grid, and *dump()* consistently stores
`header.depth = 5`, `header.n = (3,2)`. *restore()* however rebuilds
the *init_grid()* argument as a pure power of two, silently replacing
the width factor $nx$ by the next power of two, and *init_grid()*
rounds up a second time: $\log_2(2^{5+2}/3) = \log_2(42)$ rounds up to
6. The restored grid is one level deeper than the dumped one and the
traversal aborts in the multigrid `refine_cell()` stub. The two
roundings cancel exactly only when $nx$ is a power of two.

Passing *init_grid()* exactly the size the dumped grid had (the same
formula *init_grid()* itself uses to set `N`) makes
`log_base2 (N/Dimensions.x)` exact and recovers the dumped depth. This
also works with MULTIGRID_MPI (checked on 6 ranks with a 3 x 2 process
layout): there `header.n` is the process grid and `header.depth` the
local per-process depth, and the identity still holds after the
rescaling in *mpi_boundary_new()*.

~~~diff
--- a/output.h
+++ b/output.h
@@ -1265,11 +1265,7 @@
 #endif // MULTIGRID_MPI
   dimensions (header.n.x, header.n.y, header.n.z);
-  double n = header.n.x;
-  int depth = header.depth;
-  while (n > 1)
-    depth++, n /= 2;
-  init_grid (1 << depth);
+  init_grid ((1 << header.depth)*header.n.x);
 #endif // multigrid
~~~

I believe this also fixes
[my_report_bug.c](/sandbox/bugs/my_report_bug.c), which reports the
same abort with *multigrid3D* + MPI at particular process counts
(e.g. `dimensions (4,1,1)` on 256 CPUs): there the width factor seen by
*restore()* is the process grid rather than the domain aspect ratio,
but the *init_grid()* argument is reconstructed the same wrong way.

*/
