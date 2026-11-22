/**
# CUDA `DOUBLE_PRECISION=1` leaves `dvec2`/`dvec3` undefined and downcasts `cos/sin/exp/pow`

When compiling a 2D or 3D CUDA program with `DOUBLE_PRECISION=1`, any
GPU `foreach()` kernel that declares a `coord` fails NVRTC compilation
with `identifier "dvec3" is undefined` (or `dvec2` in 2D) and falls
back to CPU. In addition, `cos`, `sin`, `exp`, and `pow` are downcast
to `float(x)` even on CUDA where native FP64 transcendentals exist.

Reproduce directly via:

~~~bash
make cuda_fp64_dvec_undef.tst
~~~

Observed output on NVIDIA GeForce RTX 3070 (CUDA/NVRTC):

~~~text
cuda_fp64_dvec_undef.c:41: CUDA: error: identifier "dvec3" is undefined
cuda_fp64_dvec_undef.c:41: CUDA: !coord v = {x, y, z};
cuda_fp64_dvec_undef.c:41: CUDA: !^
cuda_fp64_dvec_undef.c:40: warning: foreach() done on CPU (see CUDA errors above)
val = 0 ref = 1.117388885112192 diff = 1.117e+00
Assertion `fabs (val - ref) < 1e-12' failed.
~~~
*/

#define DOUBLE_PRECISION 1
#include "grid/cuda/multigrid3D.h"

int main()
{
  init_grid (4);
  scalar a[];
  foreach() {
    coord v = {x, y, z};
    a[] = v.x/L0 + cos (0.123456789012345);
  }
  double ref = 0.125 + cos (0.123456789012345);
  double val = interpolate (a, 0.125, 0.125, 0.125);
  fprintf (stderr, "val = %.16g ref = %.16g diff = %.3e\n",
           val, ref, fabs (val - ref));
  assert (fabs (val - ref) < 1e-12);
}

/**
## Root cause

In [/src/grid/gpu/grid.h](/src/grid/gpu/grid.h):

1. Lines 485–491 (`#if _CUDA`) define `struct vec2` and `struct vec3`
   with `float` members, but never define `struct dvec2` or
   `struct dvec3`:

   ~~~c
   "struct vec2 { float x, y; };\n"
   "struct vec3 { float x, y, z; };\n"
   ~~~

2. Lines 520–528 (`#else // DOUBLE_PRECISION`) then alias
   `GPU_COORD2_TYPE "dvec2"`, `GPU_COORD3_TYPE "dvec3"`, and downcast
   `cos(x)`, `sin(x)`, `exp(x)`, `pow(x,y)` to `float(x)` even when
   `_CUDA` is defined (that downcast is only needed for GLSL).

## Proposed fix

1. In `/src/grid/gpu/grid.h` under `#if _CUDA`, also define:

   ~~~c
   "struct dvec2 { double x, y; };\n"
   "struct dvec3 { double x, y, z; };\n"
   ~~~

2. In `/src/grid/gpu/grid.h` under `#else // DOUBLE_PRECISION`, guard
   the `float(x)` transcendental downcasts with `#if !_CUDA` (leaving
   native FP64 `cos`, `sin`, `exp`, `pow` active on CUDA).
*/
