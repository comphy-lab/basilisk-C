/**
# OpenGL GPU reduction scratch SSBO overflows when `Dimensions.y > Dimensions.x`

On an OpenGL GPU multigrid (`-grid=gpu/multigrid`), any parallel
reduction over at least `nwgr*stride = 4096` cells aborts with
`SIGABRT` when `Dimensions.y > Dimensions.x` (or `Dimensions.z > Dimensions.x`
in 3D):

~~~text
gpu_reduct_ssbo_overflow.gpu: reduction.c:18: cpu_reduction: Assertion `a' failed.
~~~

Compile and run on GPU (aborts directly on the server):

~~~bash
make gpu_reduct_ssbo_overflow.tst
~~~

On the CPU, the reduction produces `area = 4 (expected 4)`. On the GPU
(tested on NVIDIA GeForce RTX 3070), `glMapBufferRange()` fails with
`GL_INVALID_VALUE` and returns `NULL`, triggering `assert (a)` in
`cpu_reduction()`.
*/

#include "grid/gpu/multigrid.h"
#include "utils.h"

int main()
{
  dimensions (1, 4);
  size (1. [0]);
  origin (0., 0.);
  init_grid (64); // 64 x 256 = 16384 >= 4096 cells

  double area = 0.;
  foreach (reduction(+:area))
    area += dv();
  fprintf (stderr, "area = %g (expected 4)\n", area);

  assert (fabs (area - 4.) < 1e-4);
}

/**
## Root cause

In [/src/grid/gpu/reduction.c](/src/grid/gpu/reduction.c) (lines 55–58),
the two scratch reduction SSBOs `data->reduct[0..1]` are allocated with:

~~~c
size_t size = ((dim == 2 ? sq((size_t)N + 1) : cube((size_t)N + 1))
                                            /stride + 1)*sizeof(real);
~~~

Because `N = (1 << depth())*Dimensions.x` scales only along the $x$-axis,
for `dimensions(1, 4)` and `init_grid(64)` we have $N = 64$ while the
grid has $64 \times 256 = 16384$ cells. The scratch SSBO is sized for
only $(65^2 / 64 + 1) = 67$ `real`s (`268` bytes), whereas the reduction
compute shader writes $16384 / 64 = 256$ `real`s (`1024` bytes) and
`cpu_reduction()` calls `glMapBufferRange(..., 0, 1024, GL_MAP_READ_BIT)`,
which exceeds the buffer size and returns `NULL`.

## Proposed fix

Size `data->reduct[0..1]` using `nmax = N/Dimensions.x` multiplied by
`Dimensions.{x,y,z}` along each axis:

~~~diff
--- a/src/grid/gpu/reduction.c
+++ b/src/grid/gpu/reduction.c
@@ -52,8 +52,11 @@
   if (!br[0]) {
     GL_C (glGenBuffers (2, br));
     for (int i = 0; i < 2; i++) {
       GL_C (glBindBuffer (GL_SHADER_STORAGE_BUFFER, br[i]));
-      size_t size = ((dim == 2 ? sq((size_t)N + 1) : cube((size_t)N + 1))
-                                                  /stride + 1)*sizeof(real);
+      size_t nmax = (size_t)(N/Dimensions.x);
+      size_t total = (nmax*Dimensions.x + 1)*(nmax*Dimensions.y + 1);
+      if (dim == 3)
+        total *= (nmax*Dimensions.z + 1);
+      size_t size = (total/stride + 1)*sizeof(real);
       assert (size < GPUContext.max_ssbo_size); // must fit within a single SSBO
       GL_C (glBufferData (GL_SHADER_STORAGE_BUFFER, size, NULL, GL_DYNAMIC_READ));
     }
~~~

## Related note (`gpu_init_grid`)

In [/src/grid/gpu/grid.h](/src/grid/gpu/grid.h) (lines 1744–1747),
`gpu_init_grid(int n)` starts with `if (grid && n == N) return;`,
which skips reallocating the GPU grid and SSBOs if `dimensions()` is
changed between consecutive `init_grid(n)` calls with the same `n`:

~~~diff
--- a/src/grid/gpu/grid.h
+++ b/src/grid/gpu/grid.h
@@ -1743,8 +1743,10 @@
 
 void gpu_init_grid (int n)
 {
-  if (grid && n == N)
+  static ivec prev_dim = {0};
+  if (grid && n == N && !memcmp (&Dimensions, &prev_dim, sizeof (Dimensions)))
     return;
+  prev_dim = Dimensions;
   gpu_free_grid();
   init_grid (n);
~~~
*/
