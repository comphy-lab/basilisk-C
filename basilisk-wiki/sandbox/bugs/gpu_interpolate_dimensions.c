/**
# `interpolate()` / `foreach_point()` uses `N` instead of `(nx, ny, nz)` on non-square GPU grids

On non-square / non-cubic GPU multigrids set with *dimensions()*,
*foreach_point()* (and therefore *interpolate()*) returns wrong values
or `0` when the aspect ratio is not `1:1` (2D) or `1:1:1` (3D).

Compile and run on GPU (fails directly on the server):

~~~bash
make gpu_interpolate_dimensions.tst
~~~

On the CPU, both probes return the exact cell-centered values
$s(x, y) = x + 10 y$:

~~~text
interpolate(0.8125, 0.5625) = 6.4375 (expected 6.4375)
interpolate(1.3125, 0.5625) = 6.9375 (expected 6.9375)
~~~

On the GPU (`-grid=gpu/multigrid`, tested on NVIDIA GeForce RTX 3070),
the first probe reads the wrong cell and the second reads past the end
of the level (`0`):

~~~text
interpolate(0.8125, 0.5625) = 7.1875 (expected 6.4375)
interpolate(1.3125, 0.5625) = 0 (expected 6.9375)
gpu_interpolate_dimensions.gpu.c:53: Assertion `fabs (v1 - 6.4375) < 1e-4 && fabs (v2 - 6.9375) < 1e-4' failed.
~~~
*/

#include "grid/gpu/multigrid.h"
#include "utils.h"

int main()
{
  dimensions (2, 1);
  size (2. [0]);
  origin (0., 0.);
  init_grid (16); // 16 x 8 = 128 cells (Delta = 0.125)

  scalar s[];
  foreach()
    s[] = x + 10.*y;

  double v1 = interpolate (s, 0.8125, 0.5625, linear = false); // (i=6,  j=4) -> 6.4375
  double v2 = interpolate (s, 1.3125, 0.5625, linear = false); // (i=10, j=4) -> 6.9375
  fprintf (stderr, "interpolate(0.8125, 0.5625) = %g (expected 6.4375)\n", v1);
  fprintf (stderr, "interpolate(1.3125, 0.5625) = %g (expected 6.9375)\n", v2);

  assert (fabs (v1 - 6.4375) < 1e-4 && fabs (v2 - 6.9375) < 1e-4);
}

/**
## Root cause

In [/src/grid/gpu-multigrid.h](/src/grid/gpu-multigrid.h), `val_red_(s)`
stores per-cell reduction outputs using the level strides
`NY = Nl*Dimensions.y` (and `NZ = Nl*Dimensions.z` in 3D), where
`Nl = region->level > 0 ? 1 << (region->level - 1) : N/Dimensions.x`
(matching `csOrigin` computed in the kernel dispatch).

However, in `gpu_reduction()` across all five GPU backends, the
`foreach_point` readback still bounds-checks against `N` and indexes
with stride `N` along every axis:

* [/src/grid/gpu/reduction.c](/src/grid/gpu/reduction.c) (lines 142–148)
* [/src/grid/cuda/cuda.c](/src/grid/cuda/cuda.c) (lines 756–762)
* [/src/grid/opencl/opencl.c](/src/grid/opencl/opencl.c) (lines 819–825)
* [/src/grid/hip/hip.c](/src/grid/hip/hip.c) (lines 846–852)
* [/src/grid/metal/metal.c](/src/grid/metal/metal.c) (lines 803–809)

Consequently:

* When `Dimensions.x > Dimensions.y` (wide domain, e.g. `dimensions(2, 1)`),
  `i*N + j` uses stride `N = 16` instead of `ny = 8`, reading the wrong
  cell or past the end of the level (`0`).
* When `Dimensions.y > Dimensions.x` (tall domain), points with
  $y > Y_0 + L_0$ have `j >= N`, fail `j < N`, and silently return `0`.

## Proposed fix

Compute `Nl`, `nx`, `ny`, `nz` in `gpu_reduction()` and index with
`(size_t)i*ny + j` (2D) or `((size_t)i*ny + j)*nz + k` (3D). For
[/src/grid/gpu/reduction.c](/src/grid/gpu/reduction.c):

~~~diff
--- a/src/grid/gpu/reduction.c
+++ b/src/grid/gpu/reduction.c
@@ -140,12 +140,14 @@
     GL_C (glBindBufferBase (GL_SHADER_STORAGE_BUFFER, 1 + i - start, ssbo[i]));
 
   if (is_foreach_point) {
-    real result = 0.;
-    int i = (region->p.x - X0)/L0*N;
-    int j = (region->p.y - Y0)/L0*N;
-    int k = dim == 3 ? (region->p.z - Z0)/L0*N : 0;
-    if (i >= 0 && i < N && j >= 0 && j < N && (dim == 2 || (k >= 0 && k < N))) {
-      offset += dim == 2 ? i*N + j : (i*N + j)*N + k;
+    real result = (op == 'm' ? HUGE : op == 'M' ? -HUGE : 0.);
+    int Nl = region->level > 0 ? 1 << (region->level - 1) : N/Dimensions.x;
+    int nx = Nl*Dimensions.x, ny = Nl*Dimensions.y, nz = dim == 3 ? Nl*Dimensions.z : 1;
+    int i = (region->p.x - X0)/L0*nx;
+    int j = (region->p.y - Y0)/L0*nx;
+    int k = dim == 3 ? (region->p.z - Z0)/L0*nx : 0;
+    if (i >= 0 && i < nx && j >= 0 && j < ny && (dim == 2 || (k >= 0 && k < nz))) {
+      offset += dim == 2 ? (size_t)i*ny + j : ((size_t)i*ny + j)*nz + k;
       GL_C (glUniform1ui (loffset, offset));
~~~

The same `Nl` / `(nx, ny, nz)` indexing fix applies to the `foreach_point`
block of `gpu_reduction()` in `cuda.c`, `opencl.c`, `hip.c`, and `metal.c`.
*/
