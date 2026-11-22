/**
# `DOUBLE_PRECISION=1` silently truncates GPU kernel literals to FP32

When compiling a GPU program with `-DDOUBLE_PRECISION=1`, `qcc` still
appends `f` suffixes to every floating-point literal inside generated
GPU kernels (e.g. `(1.0 + 1e-12) - 1.0` is emitted as
`(1.0f + 1e-12f) - 1.0f`), and `grid/gpu/grid.h` hardcodes `pi` and
`nodata` as `float` literals (`3.14159265359f` and `(1e30f)`).

Compile and inspect the generated kernel, or run on the GPU:

~~~bash
qcc -source gpu_fp64_literal_truncation.c
grep -n "1e-12\|3\.14159\|1e30f" _gpu_fp64_literal_truncation.c
make gpu_fp64_literal_truncation.tst
~~~

Inspecting `_gpu_fp64_literal_truncation.c` produced by `qcc -source`
shows that the GPU kernel string was emitted with `f` suffixes while
the CPU fallback kept `double` literals:

~~~c
// GPU kernel string:
val_out_(a,0,0,0) = (1.0f + 1e-12f) - 1.0f;
#define pi 3.14159265359f
#define nodata (1e30f)

// CPU fallback:
val(a,0,0,0) = (1.0 + 1e-12) - 1.0;
~~~

Running on the GPU (`make gpu_fp64_literal_truncation.tst`) fails with:

~~~text
opengl.c:692: error: OpenGl 00000502 for 'glUniform2iv (0, 1, csOrigin);'
~~~
*/

#define DOUBLE_PRECISION 1
#include "grid/gpu/multigrid.h"

int main()
{
  init_grid (4);
  scalar a[];
  foreach()
    a[] = (1.0 + 1e-12) - 1.0;

  double v = interpolate (a, 0.5, 0.5);
  fprintf (stderr, "v = %.16e (expected 1.0000000000000000e-12)\n", v);
  assert (fabs (v - 1e-12) < 1e-15);
}

/**
## Root cause

1. In [/src/include.lex](/src/include.lex) (line 328), `qcc` scans the
   raw, unpreprocessed text of included headers for:

   ~~~lex
   ^{SP}*#{SP}*define{SP}+SINGLE_PRECISION{WS}+1{SP}*$ {
     fp32 = 1;
   }
   ~~~

   Because [/src/grid/gpu-multigrid.h](/src/grid/gpu-multigrid.h)
   contains `# define SINGLE_PRECISION 1` in its `#else` branch,
   `include.lex` matches it unconditionally even when
   `-DDOUBLE_PRECISION=1` is passed, setting `.fp32 = 1` in `qcc`
   (`ast/kernels.c` then appends `f` to every floating-point literal in
   every generated GPU kernel).

2. In [/src/grid/gpu/grid.h](/src/grid/gpu/grid.h) (lines 618–619), `pi`
   and `nodata` are hardcoded with `f` suffixes regardless of
   `SINGLE_PRECISION`:

   ~~~c
   "#define pi 3.14159265359f\n"
   "#define nodata (1e30f)\n"
   ~~~

## Proposed fix

1. In `/src/include.lex` (or `/src/qcc.c`), clear `fp32 = 0` when
   `-DDOUBLE_PRECISION` is defined (or evaluate `DOUBLE_PRECISION`
   before setting `fp32`).
2. In `/src/grid/gpu/grid.h`, guard `pi` and `nodata` with
   `#if SINGLE_PRECISION ... #else ... #endif`:

   ~~~c
   #if SINGLE_PRECISION
     "#define pi 3.14159265359f\n"
     "#define nodata (1e30f)\n"
   #else
     "#define pi 3.14159265358979323846\n"
     "#define nodata (1e30)\n"
   #endif
   ~~~
*/
