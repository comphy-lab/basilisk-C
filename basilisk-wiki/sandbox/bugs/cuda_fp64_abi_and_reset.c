/**
# `libbuda.a` and `libgpu.a` hardcode FP32 (`SINGLE_PRECISION`), breaking `DOUBLE_PRECISION=1`

When compiling a CUDA or OpenGL GPU program with `DOUBLE_PRECISION=1`,
the application compiles with `sizeof(real) == 8`, while the prebuilt
runtime libraries `libbuda.a` (`grid/cuda/cuda.c`) and `libgpu.a`
(`grid/gpu/opengl.c`) are always compiled with `-DSINGLE_PRECISION` and
hardcoded `typedef float real;` (`sizeof(real) == 4`).

This produces a silent host/device ABI mismatch (`sum = 0` instead of
`384` on both CUDA and OpenGL). In addition, `cuda.c` leaves
`reset_scalar()` unimplemented for `DOUBLE_PRECISION` and hardcodes
`float` inside `cuda_reduce()`.

Reproduce directly via:

~~~bash
make cuda_fp64_abi_and_reset.tst
~~~

Observed output on NVIDIA GeForce RTX 3070 (CUDA):

~~~text
cuda_fp64_abi_and_reset.c:36: warning: CUDA: found FP64 assembly in single precision mode
sum = 0 (expected 384, diff = 3.840e+02)
cuda_fp64_abi_and_reset.c:41: Assertion `fabs (sum - expected) < 1e-12' failed.
~~~
*/

#define DOUBLE_PRECISION 1
#include "grid/cuda/multigrid.h"

int main()
{
  init_grid (16);
  scalar a[];
  reset ({a}, 1.5);

  double sum = 0.;
  foreach (reduction(+:sum))
    sum += a[];

  double expected = 1.5*sq(N);
  fprintf (stderr, "sum = %.16g (expected %.16g, diff = %.3e)\n",
           sum, expected, fabs (sum - expected));
  assert (fabs (sum - expected) < 1e-12);
}

/**
## Root cause

1. [/src/grid/cuda/Makefile](/src/grid/cuda/Makefile) (line 3) and
   [/src/grid/gpu/Makefile](/src/grid/gpu/Makefile) (line 3) hardcode
   `-DSINGLE_PRECISION`.
2. [/src/grid/cuda/cuda.c](/src/grid/cuda/cuda.c) (line 42) and
   [/src/grid/gpu/opengl.c](/src/grid/gpu/opengl.c) (line 20) hardcode
   `typedef float real;` (and `cuda.c` omits `#define HUGE 1e30f` /
   `1e30`, which also fails to compile on newer glibc).
3. [/src/grid/cuda/cuda.c](/src/grid/cuda/cuda.c) (lines 345–347)
   leaves `reset_scalar()` unimplemented under `DOUBLE_PRECISION`
   (`"error: not implemented yet"`), and `cuda_reduce()` (lines
   596–739) hardcodes `float` buffers, kernel types, and return type.

## Proposed fix

1. Guard `typedef float real;` and `HUGE` with
   `#if SINGLE_PRECISION ... #else typedef double real; #define HUGE 1e30 #endif`
   in both `cuda.c` and `opengl.c`.
2. Implement the `DOUBLE_PRECISION` branch of `reset_scalar()` and
   switch `cuda_reduce()` from `float` to `real` in `cuda.c`.
3. Build companion FP64 libraries (`libbuda64.a`, `libgpu64.a`) without
   `-DSINGLE_PRECISION` and select them via `#if DOUBLE_PRECISION` in
   `#pragma autolink` in `grid/{cuda,gpu}/multigrid{,3D}.h`.
*/
