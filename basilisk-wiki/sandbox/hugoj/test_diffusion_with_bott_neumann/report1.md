# Report: `simple.c` vertical diffusion test — CPU vs GPU, single vs double precision

*Date: 18 September 2026*

## Test setup

Source: [`simple.c`](simple.c) in this directory. A minimal layered-solver test of
`vertical_diffusion()` from [`src/layered/diffusion.h`](/src/layered/diffusion.h):

- `N = 1`, `L0 = 1`, `nl = 3` layers, each of thickness `h = 1/3`
- Tracer `T = 1` in every layer
- One call: `vertical_diffusion (point, h, T, dt=1, D=1, dst=0, sb=0, lambda_b=HUGE)`
  (Neumann top, Navier-slip bottom, both homogeneous)
- Prints `T[] - 1.` per layer; the exact answer is **0** for all layers

Build variants used:

```bash
make simple.tst                                   # CPU double (default)
CFLAGS="-DSINGLE_PRECISION" make simple.tst       # CPU float
make simple.gpu.tst                               # OpenGL float (GPU default)
make simple.cuda.tst                              # CUDA float
make simple.hip.tst                               # HIP float
make simple.ocl.tst                               # OpenCL float
CFLAGS='-DDOUBLE_PRECISION' make simple.gpu.tst  # OpenGL double
CFLAGS='-DDOUBLE_PRECISION' make simple.cuda.tst # CUDA double
```

## Results

| build | precision | errors (T−1) | in float ulps\* |
|---|---|---|---|
| CPU | double | 0.000000000000 (all layers) | 0 |
| CPU | float | 0.000000000000 (all layers) | 0 |
| OpenGL (GLSL) | float | 4.768e-7 ×3 | 4 |
| CUDA | float | 5.960e-7 ×3 | 5 |
| HIP | float | 5.960e-7 ×3 | 5 |
| OpenCL | float | 7.153e-7 ×3 | 6 |
| OpenGL (GLSL) | double | −0.000000000000 ×3 | 0 (**but ran on CPU**, see Q2) |
| CUDA | double | −1.000000000000 ×3 (**broken**, see Q2) | — |

\* 1 ulp of a float near 1.0 is 2⁻²³ ≈ 1.1920929e-7. Every single-precision error
is a small integer multiple of it.

Two questions were asked about these numbers; both are answered below.

---

## Q1: Why are single-precision results not the same across the board?

**Short answer: they are all correct to within a few float ulps; the differences
come from per-backend floating-point rounding (fast-math, FMA contraction,
instruction ordering). This is inherent to the test, not a bug.**

Details:

1. **The test is ulp-amplifying.** The exact answer is `T = 1.0`, so *any*
   single rounding difference anywhere in the ~30-operation tridiagonal
   (Thomas) solve surfaces as a small integer error in ulps: 2–3 ulp on CPU
   double-path history, 4 ulp on GLSL, 5 ulp on CUDA/HIP, 6 ulp on OpenCL.
   (Earlier runs also showed 2/3/4 ulp variations depending on build flags.)
2. **CPU** compiles with strict IEEE semantics (`gcc -O2`, no fast-math). For
   these particular values the elimination happens to cancel exactly, giving
   `1.0f` and an error of exactly 0.
3. **CUDA** kernels are compiled at runtime by NVRTC with
   `-use_fast_math` ([`src/grid/cuda/cuda.c:173`](/src/grid/cuda/cuda.c)).
   The cached PTX for this very kernel (`/tmp/buda/74d87fe0`, entry
   `init_0_53`) contains `div.approx.ftz.f32` — approximate division
   (~2 ulp each) plus denormal flushing. Hence 4–5 ulp.
4. **GLSL / OpenCL / HIP** each let the vendor driver choose its own FMA
   contraction and instruction ordering, so each lands on a slightly different
   integer number of ulps.
5. **Consequence:** identical results across backends are only achievable in
   double precision, or by testing with a tolerance. For a float test of this
   kind, compare with e.g. `fabs(T[] - 1.) < 1e-6`.

## Q2: Double precision "works" on OpenGL but not on CUDA — what do the FP64 warnings mean?

**Short answer: neither backend actually ran the test on the GPU in double.
OpenGL silently fell back to the CPU (due to a typo bug in the GLSL
preprocessor definitions) and thus printed the CPU's exact zero. CUDA ran a
double kernel against a single-precision-only runtime library, which corrupts
the host↔device transfer of `T` (it comes back as 0, hence `T−1 = −1`). The
FP64 warning is precisely the detector of that CUDA mismatch.**

### OpenGL double = silent CPU fallback

The double run prints:

```
(fragment shader):482: GLSL: error C7101: Macro coord redefined
simple.gpu.c:17: warning: foreach() done on CPU (see GLSL errors above)
(fragment shader):482: GLSL: error C7101: Macro coord redefined
simple.gpu.c:23: warning: foreach() done on CPU (see GLSL errors above)
```

In double precision and `dimension == 2`, the GLSL preprocessor string in
[`src/grid/gpu/grid.h`](/src/grid/gpu/grid.h) emits **two** definitions of
`coord`:

- line 477: `#define coord dvec3` (general double case)
- line 487: `#define coord dvec2` (the `dimension == 2`, `!SINGLE_PRECISION`
  branch — compare line 485 in the single-precision branch, which correctly
  defines **`_coord`**)

Line 487 should read `"#define _coord dvec2\n"` — it lost its underscore.
GLSL refuses the macro redefinition, shader compilation fails, and Basilisk
falls back to running both `foreach()` loops on the **CPU** (in double), which
prints the exact `−0.000000000000`. The GPU never computed anything. This bug
is present in the upstream (darcs) version of `grid.h`, unmodified locally.

### CUDA double = double kernel vs single-precision runtime library

The CUDA driver library `libbuda.a` is built with `-DSINGLE_PRECISION`
hardcoded ([`src/grid/cuda/Makefile:3`](/src/grid/cuda/Makefile)), and
`cuda.c` hardcodes `typedef float real` (line 42) regardless of
`DOUBLE_PRECISION`. When the *main program* is built with
`-DDOUBLE_PRECISION`:

- the runtime-generated kernel source says `#define real double`
  (from `gpu/grid.h:476`), so NVRTC compiles **double** kernels (`.f64` PTX);
- but all host↔device sizing in libbuda uses `sizeof(real)` = 4 bytes:
  - `gpu_cpu_sync_scalar` (cuda.c:325–327) copies the wrong sizes/offsets
    for 8-byte host fields → `T` is never correctly transferred back;
  - `reset_scalar` (cuda.c:338) has a double branch that is literally
    unimplemented: `fprintf (stderr, "%s:%d: error: not implemented yet\n")`
    (cuda.c:347–350) — and that `fprintf` is even missing its arguments;
- result: the kernel runs in double on the device, but the host reads back
  zeros → `T−1 = −1`.

### The FP64 warning, decoded

Emitted by `compile_ptx` under `#if SINGLE_PRECISION`
([`src/grid/cuda/cuda.c:202–207`](/src/grid/cuda/cuda.c)) when the compiled
PTX contains `.f64` instructions:

> "the runtime-generated kernel is double, but this buda library was built
> single-precision"

So the warnings at `simple.c:17` and `:23` were a **correct diagnosis** of
exactly the mismatch above. They are not about the source file's precision per
se — they compare the *kernel's* PTX against the *library's* build flag.

### Same structural issue elsewhere

All GPU backends hardcode single precision in their library builds:

- `src/grid/gpu/Makefile:3` (libgpu.a, OpenGL) and `opengl.c:19`
  (`typedef float real`)
- `src/grid/hip/Makefile:3`
- `src/grid/opencl/Makefile:3`

So GPU double is only partially implemented in this tree: the **main
program** honors `DOUBLE_PRECISION` (via `gpu-multigrid.h` /
`gpu-cartesian.h`, which map it to `SINGLE_PRECISION=0`), but the
**precompiled backends do not**.

---

## Proposed fixes

1. **GLSL double (5 min):** fix the `_coord` typo at
   `src/grid/gpu/grid.h:487` (`#define coord dvec2` → `#define _coord dvec2`).
   Double shaders then compile and the CPU fallback disappears.
2. **CUDA double (2–4 h):** make `cuda.c` precision-consistent:
   - build `libbuda.a` without the hardcoded `-DSINGLE_PRECISION`
     (derive it from `DOUBLE_PRECISION` like the main program does);
   - guard `typedef float real` → `double` when `DOUBLE_PRECISION`;
   - implement the `reset_scalar` double branch (`cuMemsetD64`), fixing the
     malformed `fprintf`.
   With consistency, the FP64 warning logic also becomes meaningful again.
3. **OpenGL double (1–2 h):** same treatment for `libgpu.a` (`opengl.c`
   `typedef` + its `#if SINGLE_PRECISION` uniform-upload blocks; the double
   branches already exist, e.g. `glUniform1dv` at opengl.c:655–662).
4. **HIP / OpenCL (1–2 h each):** same recipe afterwards.
5. **Verification:** `CFLAGS='-DDOUBLE_PRECISION' make simple.cuda.tst` and
   `simple.gpu.tst` should print ~1e-16 (or exactly 0) with no FP64 warning;
   keep single-precision runs as a tolerance test (< 1e-6), per Q1.

## Side notes

- `src/layered/diffusion.h` currently contains an uncommitted local edit
  (commented out): `//@define double float` around `vertical_diffusion()`
  (lines 31/125). It is unnecessary (precision is controlled by
  `SINGLE_PRECISION`/`DOUBLE_PRECISION` which switch the `real` typedef) and
  ineffective on GPU backends — recommend removing it.
- The upstream test `src/test/layers.cuda` fails in this tree independently of
  the above (2D/`LAYERS` parse errors in the GPU kernel generation,
  `Point.l`/`block` handling). Worth tracking separately before relying on
  layered GPU runs.
- Local edits force `--gpu-architecture=compute_89` / `sm_89`
  (cuda.c:148, hip.c:181–184) for a GPU newer than the supported arch list.
  Harmless here, but they are debug leftovers.





  Follow up:
 Short answer: No — the errors do not come from the top layer. Your new data actually proves the opposite: for h=1 and h=2 the error appears at layer 1 (near the bottom), while for h=1/3 it appeared at layers 4–5. If the top-BC row were the culprit, the error would stay at nl−1 regardless of h. The position moves with h because it comes from rounding in the Thomas-algorithm divisions (diffusion.h lines 116–123), not from the boundary rows.
What I verified in diffusion.h (nl=10, h=1)
With dt=1, D=1, dst=0, s_b=0, lambda_b=HUGE=1e30, every coefficient is exact in fp32:
- interior (57–61): a=c=−1, b=3
- top (80–82): a₉=−1, b₉=2
- bottom (100–105): den = 1·4 + 6e30 → b₀=2.0, and the tiny correction terms (7/den ~ 6e-31) vanish in float → c₀=−1
So the exact system is [[2,−1],[−1,3,−1],…,[−1,3,−1],[−1,2]] · s = 1 → s=1 identically. The only inexact operations left are the divisions in
b[l]   -= a[l]*c[l-1]/b[l-1];    // line 117
rhs[l] -= a[l]*rhs[l-1]/b[l-1];  // line 118
s = rhs/b                          // lines 120, 123
In exact arithmetic the eliminated diagonals form a Fibonacci-ratio chain b_l → φ² ≈ 2.618 — not float-representable, so each division rounds.
Why your numbers look the way they do
- CPU fp32 (all zeros): IEEE divisions are exactly rounded; for this matrix the small perturbations from each division cancel (tie-to-even), giving exactly 1.0f at every layer.
- CUDA (±½ ulp at layer 1): kernels are compiled with -use_fast_math (cuda.c:173), so divisions become div.approx.ftz.f32 — approximate, up to ~2 ulp each. I confirmed this instruction in the cached PTX for this kernel. One of the ~30 divisions lands 1 ulp off at layer 1, then the back-substitution contracts it by |c/b|≈0.38 per layer (and s₁ at 1−2⁻²⁴ ties to even → s₀ recovers exactly 1.0).
- h=1/3 errors at 4–5 vs h=1,h=2 error at 1: each h gives a different exact system (a=c=−D·dt/h), hence a different float error-chain → the first-failure layer moves. Not related to nl−1.
Note on the two outputs you pasted
Both blocks (labeled h=1 and h=2) are byte-identical. Either the run really gives the same result for both, or the second is a copy of the first — worth a quick re-check that the binary was rebuilt between the two (h=2 should re-trigger simple.s).
From earlier check, real quirks in diffusion.h (unrelated to this ulp):
1. nl == 1 reads h[0,0,-1] (out-of-bounds, line 80) — unused result, but UB.
2. c[nl-1] never initialized (line 35) — mathematically unused, fragile.
3. nl == 1 double-adds dst then cancels it (107–110) — convoluted.
Plan (on your approval)
1. Pinpoint the exact division: build a debug variant printing b[l], rhs[l], s[l] as hex floats on CPU + CUDA (~15 min).
2. Test the hypothesis that fast-math is the cause: recompile the CUDA backend without -use_fast_math (temporary edit of cuda.c opts) and see if CUDA then gives all zeros like CPU. (~20 min; informs whether you'd want a "safe-math" build mode.)
3. Keep the runtest comparison tolerance-based (< 2e-6 covers all observed float-backends).
