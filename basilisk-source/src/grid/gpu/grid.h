/**
# Grids on GPUs

The files in this directory implement Cartesian and Multigrid grids on
[Graphics Processing
Units](https://en.wikipedia.org/wiki/Graphics_processing_unit)
(GPUs). The ultimate goal is to allow running any [Basilisk
solver](/src/README#solvers) on GPUs without any modification to the
original source code.

To do so the [Basilisk preprocessor](/src/ast/README) automatically
generates "computation kernels" for each [loop
iterator](/Basilisk%20C#iterators). These kernels are then dynamically
compiled (at runtime) by the [OpenGL Shading
Language](https://en.wikipedia.org/wiki/OpenGL_Shading_Language)
(GLSL) compiler which is part of the (OpenGL) graphics card driver. If
compilation is successful, the corresponding loop is performed on the
GPU, otherwise the CPU is used. If this hybrid GPU/CPU hybrid mode of
operation is used, synchronisation between the GPU and CPU memory is
necessary and is done automatically.

OpenGL is an open standard (unlike
e.g. [CUDA](https://en.wikipedia.org/wiki/CUDA)) and is widely
supported by graphics cards (with the notable exception of Apple
graphics cards and some high-end "professional" Nvidia cards).

## Running on GPUs

As described above, from a Basilisk perspective GPUs are just another
type of grid. Selecting a "GPU grid" can simply be done using either

~~~literatec
#include "grid/gpu/multigrid.h"
~~~

in the source code, or using the `-grid` command line option of
[qcc](/src/qcc.c) like this

~~~bash
qcc -autolink -Wall -O2 -grid=gpu/multigrid code.c -o code -lm
~~~

The standard Basilisk [Makefile](/Tutorial#using-makefiles) also
includes the handy recipe

~~~bash
make code.gpu.tst
~~~

which will compile and run `code.c` using the `gpu/multigrid` grid.

Note that for all this to work properly you first need to
[install](#installation) the Basilisk GPU libraries.

## Installation

Basilisk uses the [GLFW](https://www.glfw.org/) library to configure
and access the graphics card and [OpenGL](https://www.opengl.org/)
(version >= 4.3) for the rest. These libraries and the associated
Basilisk libraries can be easily installed on Debian-like systems
using

~~~bash
sudo apt install libglfw3-dev
cd $BASILISK/grid/gpu
make
~~~

Note that you will also need the appropriate graphics card drivers
(often proprietary for Nvidia). Note also that (reasonably high-end)
laptop computers often have two graphics cards: a low-power, slow one
and a high-power, fast one. To check which one you are currently using
you can use something like

~~~bash
sudo apt install mesa-utils
glxinfo -B
~~~

On my Dell XPS laptop I can switch to the (proprietary driver of the)
fast Nvidia graphics card using

~~~bash
__NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia glxinfo -B
~~~

## Tests

There are several [test cases for GPUs](/src/test/README#gpu-tests)
you can try. For example

~~~bash
cd $BASILISK/test
CFLAGS=-DPRINTNSHADERS make gpu.gpu.tst
~~~

If this worked, you can then try a more interesting example

~~~bash
CFLAGS='-DSHOW' make bump2D-gpu.tst
__NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia ./bump2D-gpu/bump2D-gpu 10
~~~

and also

~~~bash
cd $BASILISK/examples
CFLAGS='-DSHOW -DBENCHMARK' make turbulence.gpu.tst
__NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia ./turbulence.gpu/turbulence.gpu 1024
~~~

## Writing code compatible with GPUs

GPUs are fast compared to CPUs because they use specialised hardware
which relies on highly-parallel (tens of thousands of execution
threads) asynchronous accesses to fast video memory channels. This
imposes strong constraints on programs which can run efficiently on
these systems, in particular regarding memory allocation and
accesses. These constraints are reflected in the programming languages
usable on GPUs, for example the [OpenGL Shading
Language](https://www.khronos.org/opengl/wiki/Core_Language_(GLSL))
(GLSL) which underlies the GPU grid in Basilisk.

GLSL is mostly a subset of C99 and the good news is that this subset
happens to be what is used within most `foreach` loops in Basilisk
(this is not a coincidence...). Thus, in many cases, simple and
efficient Basilisk code will also run transparently and efficiently on
GPUs.

### GPU/CPU hybrid mode

There are obvious cases where foreach loops will not run on GPUs (see
the [next section](#what-cannot-be-done-on-gpus)). In theses cases,
Basilisk will automatically switch to running the loop on the CPU and
will synchronize the CPU and GPU memories. Note that this also means
that the memory (i.e. scalar fields etc.) for a program is always
allocated twice: once on the CPU and once on the GPU.

As an example, consider the following simple code

~~~literatec
int main()
{
  init_grid (16);
  scalar s[];
  double k = 2.;
  foreach()
    s[] = cos(k*x)*sin(k*y);
  foreach()
    printf ("%g %g %g\n", x, y, s[]);
}
~~~

this can be run on the CPU using e.g.

~~~bash
make test.tst
~~~

If we now run on the GPU using

~~~bash
make test.gpu.tst
~~~

we get

~~~bash
test.gpu.c:9: GLSL: error: unknown function 'printf'
test.gpu.c:8: warning: foreach() done on CPU (see GLSL errors above)
~~~

Basilisk warns us that "printf" is not known in GLSL (at line 9)
and that, as a consequence, the loop at line 8 (i.e. the second loop
which includes "printf") was run on the CPU. Note that the first
message is a "GLSL: error" but that the code still ran OK on the
CPU. Note also that this error happened at runtime and not during
compilation. That's because foreach loops are compiled dynamically at
runtime by the graphics card driver.

Since GPUs have a very limited access to the operating system
(i.e. only through the OpenGL interface) we cannot expect the loop
including "printf" (or any other output) to run on the GPU. Note also
that the second loop should be "serial" rather than parallel (see
[Parallel Programming](/Basilisk C#parallel-programming)). So we need
to modify the code to

~~~literatec
...
  foreach (serial)
    printf ("%g %g %g\n", x, y, s[]);
...
~~~

If we now recompile and run with `make test.gpu.tst`, the GLSL error
and warnings are gone since we explicitely specified that the second
loop should run on the CPU (and in serial mode).

Another way to specify that a given loop should run on the CPU (either
in serial or parallel mode) is to use

~~~literatec
foreach (cpu)
  ...
~~~

Similarly one could use `foreach (gpu)` to force running on the
GPU, in which case the GLSL warning above would become an error. This
can be useful when debugging GPU codes and used in combination with
the `-cpu` [compilation flag](/src/qcc.c) which will force loops to
run on the CPU by default.

### Variable-size arrays

In C99 variable-size arrays can be defined simply using for example

~~~c
void func1 (int n, double a[n]) {
  ...
}
...
{
  int m = ...;
  double b[m];
  func (m, b);
}
~~~

Since this relies on dynamic memory allocation on the stack, this is
not possible in general in GLSL. The only cases where this will work
is if the size of the array can be computed "statically" i.e. at the
time the GLSL kernel is compiled. Furthermore, the GLSL compiler is
strict (or not very clever) and a code looking like

~~~c
int n = 3;
double a[n];
~~~

will fail with an error like

~~~bash
GLSL: error: array size must be a constant valued expression
~~~

To fix this one needs to write instead

~~~c
const int n = 3;
double a[n];
~~~

Note that the size of the array must be a constant, but only at the
time when the GLSL kernel is compiled. This allows using variable-size
arrays also in GLSL, provided their size is constant within the
kernel. For example the following code will work fine on the GPU, even
if `n` changes between calls to `func2()`.

~~~literatec
void func2 (const int n) {
  foreach() {
    const int size = n + 1;
    double a[size];
    ...
  }
}
~~~

Finally, using variable-sized arrays as function parameters, as
done in `func1()` above, is not allowed in GLSL. To work around this
strong limitation, the [kernel preprocessor](/src/ast/kernels.c) will
expand calls to functions using variable-size arrays (using the
[macro](/src/ast/macro.h) engine). Note that this means that the
function must respect the constraints applying to macros, in
particular they [can return only at the end of the
function](/src/ast/macro.h#complex-return-macros).

### What cannot be done on GPUs

* Inputs/Outputs: The only possible direct output on GPUs is the screen (see
  [output_ppm on GPUs](output.h)). All other inputs or outputs must go through
  the CPU.
* Complex memory allocation and access: There is no notion of "memory
  stack" on GPUs, all memory requirements are ***static*** and must be
  defined at compilation time. This means that variable/dynamical
  arrays, dynamic memory allocation (malloc etc.) and ***pointers***
  do not exist on GPUs (and in GLSL). Using any of these in a foreach
  loop will give you a GLSL error as above.
* Limited support for function pointers: function pointers are
  fundamentally different from memory pointers.  Basilisk includes
  limited support for function pointers i.e. what is sufficient for
  *field attributes* implementing custom boundary conditions or
  coarsening/refinement functions.
* Using external libraries: GPUs cannot (obviously) use functions defined in
  external (CPU) libraries.

## Current limitations

Aside from the fundamental constraints above, the current
implementation also has the following limitations, some of which will
be lifted in the (not-too-distant) future. In rough order of "lifting
priority" these are:

* Only 2D Cartesian and Multigrid grids for now: 3D multigrid will
  follow easily, quadtrees and octrees are more difficult.
* The maximum size of any scalar field is limited to what can be
  indexed using a 32-bits unsigned integer i.e. 2^32^ floats or 16 GB.
* Boundary conditions have only been implemented for 3x3 stencils.
* At this stage only a few solvers [have been
  tested](/src/test/README#gpu-tests). Other solvers may or may not
  work. In particular [surface tension](/src/tension.h) will not work
  yet because the [estimation of curvature](/src/curvature.h) relies
  on [code which is not portable to GPUs](/src/parabola.h).
* Only simple external types (`int`, `bool`, `float`, `double`,
  `coord` etc.) are currently supported: for example custom `typedefs`
  are not supported yet for external variables. Loop-local variables
  and functions can use (locally-defined) custom types.
* Loop-local [lists](/Basilisk%20C#lists) of scalars have very limited
  support (but are not used much anyway): external loops support is
  OK.
* Double precision (64-bits floats) is supported by Basilisk (use
  `CFLAGS='-DDOUBLE_PRECISION'`) but depends on the (often limited)
  support by graphics cards and their drivers. Note also that using
  single precision can have an important impact on the convergence and
  accuracy of [multigrid solvers](/src/poisson.h).

## Performance

To convince yourself that GPUs are worth the trouble, see the [GPU
Benchmarks](Benchmarks.md): speedups of one to two orders of magnitude
compared to CPUs are achievable.

To maximize performance, here are a few tips and observations:

* Make sure that you are using the correct graphics card and driver
 (see [glxinfo](#installation) above).
* GPUs are highly parallel so will only provide speedups for large
  enough simulations (e.g. larger than 128^2^), increasingly so as
  resolution is increased.
* Frequent CPU/GPU memory synchronisation will kill performance. Be
  careful to control how often you output data for example, much more
  so than when running on CPUs. An exception is [graphical
  outputs](output.h) which are much cheaper on GPUs and can be done
  often with little overhead. Loops done on the CPU within
  e.g. timestep iterations will generally kill performance.
* Use [built-in profiling](/src/README.trace) to check where time is
  spent. Use the `-DTRACE=3` compilation flag to get profiling
  information at the level of foreach loops.

## Profiling with Nsight Compute or Nsight Systems

~~~bash
ncu --set default ./my_app

nsys profile --stats=true <your_executable> [args]
nsys-ui report*.nsys-rep
~~~
  
## Bugs

* Trying to use more than one shader storage buffer (SSBO) with Intel
  drivers does not seem work. This limits the amount of available
  video memory to a single SSBO with a maximum size which is usually
  2GB (or less).

## See also

* [Test cases on GPUs](/src/test/README#gpu-tests)
* [GPU benchmarks](Benchmarks.md)
* [Computation kernels](/src/ast/kernels.c)
* [An experimental CUDA driver](/src/grid/cuda/cuda.c)
* [An experimental HIP driver](/src/grid/cuda/hip.c)

# Implementation */

bool on_cpu = false;

#define gpu_grid ((GridGPU *)grid)

GPUContext_t GPUContext = {
  .current_shader = -1,
};

#if _CUDA
#    define DEFINITIONS              \
  "struct ivec2 { int x, y; };\n" \
  "typedef unsigned int uint;\n" \
  "struct vec2 { float x, y; };\n" \
  "struct vec3 { float x, y, z; };\n" \
  "__host__ __device__ inline int clamp(int x, int a, int b) { return max(a, min(x, b)); }\n" \
  "#define forin(type,s,list) for (int _i = 0; _i < sizeof(list)/sizeof(type) - 1; _i++) { type s = list[_i];\n" \
  "#define forin2(a,b,c,d) for (int _i = 0; _i < sizeof(c)/sizeof(a) - 1; _i++)" \
  "  { a = c[_i]; b = d[_i];\n"                                         \
  "#define forin3(a,b,e,c,d,f) for (int _i = 0; _i < sizeof(c)/sizeof(a) - 1; _i++)" \
  "  { a = c[_i]; b = d[_i]; e = f[_i];\n"                              \
  "#define layout(x)\n"
#else // !_CUDA
#  define DEFINITIONS                                                   \
  "#define forin(type,s,list) for (int _i = 0; _i < list.length() - 1; _i++) { type s = list[_i];\n" \
  "#define forin2(a,b,c,d) for (int _i = 0; _i < c.length() - 1; _i++)" \
  "  { a = c[_i]; b = d[_i];\n"                                         \
  "#define forin3(a,b,e,c,d,f) for (int _i = 0; _i < c.length() - 1; _i++)" \
  "  { a = c[_i]; b = d[_i]; e = f[_i];\n"
#endif // !_CUDA

const char glsl_preproc[] =
  "// #line " xstr(LINENO) " " __FILE__ "\n"
#if _CUDA
#if _HIP_AMD
  "#define uniform __device__\n"
#else
  "#define uniform __constant__\n"
#endif
#endif
  DEFINITIONS
#if _CUDA
  "#define neighborp(_i,_j,_k) Point{point.i+_i,point.j+_j,point.level,point.n"
#if LAYERS
  ",point.l"
#endif
  "}\n"
#else // !CUDA
  "#define neighborp(_i,_j,_k) Point(point.i+_i,point.j+_j,point.level,point.n"
#if LAYERS
  ",point.l"
#endif
  ")\n"
#endif // !_CUDA
  "#define dimensional(x)\n"
  "#define fmin(a,b) min(a,b)\n"
  "#define fmax(a,b) max(a,b)\n"
#if !SINGLE_PRECISION
  "#define real double\n"
  "#define coord dvec3\n"
#else // !SINGLE_PRECISION
  "#define real float\n"
  "#define coord vec3\n"
#endif // !SINGLE_PRECISION
  "#define ivec ivec2\n"
  "struct scalar { int i, index; };\n"
#if dimension == 2
#if SINGLE_PRECISION
  "#define _coord vec2\n"
#else
  "#define coord dvec2\n"
  "#define cos(x) cos(float(x))\n"
  "#define sin(x) sin(float(x))\n"
  "#define exp(x) exp(float(x))\n"
  "#define pow(x,y) pow(float(x), float(y))\n"
#endif
  "struct vector { scalar x, y; };\n"
  "struct tensor { vector x, y; };\n"
#endif // dimension == 2
  "#define GHOSTS " xstr(GHOSTS) "\n"
  "struct Point { int i, j, level;"
#if MULTIGRID
  " ivec2 n;"
#else
  " uint n;"
#endif
#if LAYERS
  " int l;"
#endif
  "};\n"
  "#define field_size() _field_size\n"
  "#define ast_pointer(x) (x)\n"
  GPU_CODE()
  "#define _NVARMAX 65536\n"
  "#define is_constant(v) ((v).i >= _NVARMAX)\n"
  "#define NULL 0\n"
  "#define val(s,k,l,m) ((s).i < _NVARMAX ? valt(s,k,l,m)"
  " : _constant[clamp((s).i -_NVARMAX,0,_nconst-1)])\n"
  "#define val_out_(s,i,j,k) valt(s,i,j,k)\n"
  "#define diagonalize(a)\n"
  "#define val_diagonal(s,i,j,k) real((i) == 0 && (j) == 0 && (k) == 0)\n"
  "#define _attr(s,member) (_attr[(s).index].member)\n"
  "#define endforin() }\n"
#if LAYERS
  "#define _index(a,m) ((a).i + (point.l + _layer + (m) < _attr(a,block) ? point.l + _layer + (m) : 0))\n"
#else
  "#define _index(a,m) ((a).i)\n"
#endif
  "#define endforin2() }\n"
  "#define endforin3() }\n"
  "#define NOT_UNUSED(x)\n"
  "#define pi 3.14159265359f\n"
  "#define nodata (1e30f)\n"
  "#define fabs(x) abs(x)\n"
  "const real z = 0.;\n"
  "const int ig = 0, jg = 0;\n"
  "layout (location = 0) uniform ivec2 csOrigin = {0,0};\n"
  "layout (location = 1) uniform vec2 vsOrigin = {0.f,0.f};\n"
  "layout (location = 2) uniform vec2 vsScale = {1.f,1.f};\n"
  ;

static inline int list_size (const External * i)
{
  int size = 0;
  if (i->type == sym_SCALAR) {
    size = 1;
    if (i->nd == 1)
      for (scalar s in i->pointer) size++;
  }
  else if (i->type == sym_VECTOR) {
    size = 1;
    if (i->nd == 1)
      for (vector v in i->pointer) size++;
  }
  else if (i->type == sym_TENSOR) {
    size = 1;
    if (i->nd == 1)
      for (tensor t in i->pointer) size++;
  }
  else
    return 0;
  return size;
}

static char * write_scalar (char * fs, scalar s)
{
  char i[20], index[20];
  snprintf (i, 19, "%d", s.i);
  if (s.i < 0 || is_constant(s))
    strcpy (index, "0");
  else {
    if (s.block < 0)
      s.i += s.block;
    snprintf (index, 19, "%d", s.gpu.index - 1);
  }
  return str_append (fs, "{", i, ",", index, "}");
}

static char * write_vector (char * fs, vector v)
{
  fs = str_append (fs, "{");
  fs = write_scalar (fs, v.x);
  fs = str_append (fs, ",");
  fs = write_scalar (fs, v.y);
  fs = str_append (fs, "}");
  return fs;
}

static char * write_tensor (char * fs, tensor t)
{
  fs = str_append (fs, "{");
  fs = write_vector (fs, t.x);
  fs = str_append (fs, ",");
  fs = write_vector (fs, t.y);
  fs = str_append (fs, "}");
  return fs;
}

static scalar * apply_bc_list;

static int bc_period_x = -1, bc_period_y = -1;

static void boundary_top (Point point, int i)
{
  bool data = false;
  for (scalar s in apply_bc_list)
    if (!s.face || s.i != s.v.y.i) {
      scalar b = (s.v.x.i < 0 ? s : s.i == s.v.y.i ? s.v.x : s.v.y);
      foreach_blockf(s)
	s[i,-bc_period_y] = b.boundary_top (neighborp(i), neighborp(i,-bc_period_y), s, &data);
    }
}

static void boundary_bottom (Point point, int i)
{
  bool data = false;
  for (scalar s in apply_bc_list)
    if (!s.face || s.i != s.v.y.i) {
      scalar b = (s.v.x.i < 0 ? s : s.i == s.v.y.i ? s.v.x : s.v.y);
      foreach_blockf(s)
	s[i,bc_period_y] = b.boundary_bottom (neighborp(i), neighborp(i,bc_period_y), s, &data);
    }
}

static
void apply_bc (Point point)
{
  bool data = false;
  // face BCs
  if (point.i == GHOSTS)
    for (scalar s in apply_bc_list)
      if (s.face && s.i == s.v.x.i && s.boundary_left)
	foreach_blockf(s)
	  s[] = s.boundary_left (point, neighborp(bc_period_x), s, &data);
  if (point.i == N*Dimensions.x + GHOSTS)
    for (scalar s in apply_bc_list)
      if (s.face && s.i == s.v.x.i && s.boundary_right)
	foreach_blockf(s)
	  s[] = s.boundary_right (neighborp(bc_period_x), point, s, &data);
  if (point.j == GHOSTS)
    for (scalar s in apply_bc_list)
      if (s.face && s.i == s.v.y.i) {
	scalar b = s.v.x;
	if (b.boundary_bottom)
	  foreach_blockf(s)
	    s[] = b.boundary_bottom (point, neighborp(0,bc_period_y), s, &data);
      }
  if (point.j == N*Dimensions.y + GHOSTS)
    for (scalar s in apply_bc_list)
      if (s.face && s.i == s.v.y.i) {
	scalar b = s.v.x;
	if (b.boundary_top)
	  foreach_blockf(s)
	    s[] = b.boundary_top (neighborp(0,bc_period_y), point, s, &data);
      }
  // centered BCs
  if (point.i == GHOSTS) { // left
    for (scalar s in apply_bc_list)
      if (!s.face || s.i != s.v.x.i)
	foreach_blockf(s)
	  s[bc_period_x] = s.boundary_left (point, neighborp(bc_period_x), s, &data);
    if (point.j == GHOSTS)
      boundary_bottom (point, bc_period_x); // bottom-left
    if (point.j == N*Dimensions.y + GHOSTS - 1)
      boundary_top (point, bc_period_x);    // top-left
  }
  if (point.i == N*Dimensions.x + GHOSTS - 1) { // right
    for (scalar s in apply_bc_list)
      if (!s.face || s.i != s.v.x.i)
	foreach_blockf(s)
	  s[- bc_period_x] = s.boundary_right (point, neighborp(- bc_period_x), s, &data);
    if (point.j == GHOSTS)
      boundary_bottom (point, - bc_period_x); // bottom-right
    if (point.j == N*Dimensions.y + GHOSTS - 1)
      boundary_top (point, - bc_period_x);    // top-right
  }
  if (point.j == GHOSTS)
    boundary_bottom (point, 0);  // bottom
  if (point.j == N*Dimensions.y + GHOSTS - 1)
    boundary_top (point, 0);     // top
}

static bool is_boundary_attribute (const External * g)
{
  return (g->name[0] == '.' &&
	  (!strcmp (g->name, ".boundary_left") ||
	   !strcmp (g->name, ".boundary_right") ||
	   !strcmp (g->name, ".boundary_bottom") ||
	   !strcmp (g->name, ".boundary_top")));
}

static
void hash_external (Adler32Hash * hash, const External * g, const ForeachData * loop, int indent)
{
  if (g->type == sym_function_declaration || g->type == sym_function_definition) {
    bool boundary = is_boundary_attribute (g);
    for (scalar s in baseblock)
      if (g->name[0] != '.' || (!boundary && s.gpu.index) ||
          (boundary && (s.stencil.io & s_output))) {
	void * ptr = g->name[0] != '.' ? g->pointer :
	  *((void **)(((char *) &_attribute[s.i]) + g->nd));
	if (ptr) {
	  External * p = _get_function ((long) ptr);
	  if (p && !p->used) {
	    p->used = true;
	    a32_hash_add (hash, &ptr, sizeof (void *));
	    for (const External * e = p->externals; e && e->name; e++)
	      hash_external (hash, e, loop, indent + 2);
	  }
	}
	if (g->name[0] != '.')
	    break;
      }
  }
  else if (g->name[0] == '.') {
    int size = 0;
    switch (g->type) {
    case sym_INT: size = sizeof (int); break;
    case sym_BOOL: size = sizeof (bool); break;
    case sym_FLOAT: size = sizeof (float); break;
    case sym_DOUBLE: size = sizeof (double); break;
    case sym_IVEC: size = sizeof (ivec); break;
    case sym_SCALAR: size = sizeof (scalar); break;
    case sym_VECTOR: size = sizeof (vector); break;
    case sym_TENSOR: size = sizeof (tensor); break;
    default: return;
    }
    for (scalar s in baseblock)
      if (s.gpu.index) {
	char * data = (char *) &_attribute[s.i];
	a32_hash_add (hash, data + g->nd, size);
      }
  }
  else if (g->type == sym_SCALAR || g->type == sym_VECTOR || g->type == sym_TENSOR) {
    assert (g->nd == 0 || g->nd == 1);
    void * pointer = g->pointer;
    int size;
    if (g->nd == 1) {
      size = 0;
      for (scalar s in pointer)
	size += sizeof (scalar);
    }
    else if (g->type == sym_SCALAR)
      size = sizeof (scalar);
    else if (g->type == sym_VECTOR)
      size = sizeof (vector);
    else // sym_TENSOR
      size = sizeof (tensor);
    a32_hash_add (hash, pointer, size);
  }
  else if (IS_EXTERNAL_CONSTANT(g))
    a32_hash_add (hash, g->pointer, sizeof(int));
}

static
uint32_t hash_shader (const External * externals,
		      const ForeachData * loop,
		      const RegionParameters * region, const char * kernel)
{
  Adler32Hash hash;
  a32_hash_init (&hash);
  a32_hash_add (&hash, &region->level, sizeof (region->level));
  a32_hash_add (&hash, &N, sizeof (N));
#if LAYERS
  a32_hash_add (&hash, &nl, sizeof (nl));
#endif
  a32_hash_add (&hash, &Dimensions, sizeof (Dimensions));
  a32_hash_add (&hash, &Period, sizeof (Period));
  a32_hash_add (&hash, kernel, strlen (kernel));
  a32_hash_add (&hash, &nconst, sizeof (nconst));
  a32_hash_add (&hash, _constant, sizeof (*_constant)*nconst);
  static External ext[] = {
    { .name = ".boundary_left",   .type = sym_function_declaration,
      .nd = attroffset (boundary_left) },
    { .name = ".boundary_right",  .type = sym_function_declaration,
      .nd = attroffset (boundary_right) },
    { .name = ".boundary_top",    .type = sym_function_declaration,
      .nd = attroffset (boundary_top) },
    { .name = ".boundary_bottom", .type = sym_function_declaration,
      .nd = attroffset (boundary_bottom) },
#if LAYERS
    { .name = ".block", .type = sym_INT, .nd = attroffset (block) },
#endif
    { .name = NULL }
  };
  foreach_function (f, f->used = false);
  for (const External * g = loop->dirty ? ext : ext + 4; g->name; g++)
    hash_external (&hash, g, loop, 2);
  int imax = 1;
  if (GPUContext.nssbo > 1)
    for (scalar s in baseblock)
      if (s.gpu.index && s.i + s.block > imax)
	imax = s.i + s.block;
  for (const External * g = externals; g && g->name; g++) {
    if (g->reduct) {
      a32_hash_add (&hash, &g->s, sizeof (scalar));
      if (GPUContext.nssbo > 1) {
	scalar s = g->s;
	if (s.i + s.block > imax)
	  imax = s.i + s.block;
      }
    }
    hash_external (&hash, g, loop, 2);
  }
  a32_hash_add (&hash, &imax, sizeof (imax));
  return a32_hash (&hash);
}

static
bool is_void_function (char * code)
{
  while (*code != '/') code++; code++;
  while (*code != '/') code++; code++;
  while (*code != '\n') code++; code++;
  while (strchr ("\n\r\t ", *code)) code++;
  return !strncmp (code, "void ", 5);
}

static char * type_string (const External * g)
{
  switch (g->type) {
  case sym_DOUBLE:
#if !SINGLE_PRECISION
    return "double"; break;
#endif
  case sym_FLOAT: return "float";
  case sym_function_declaration: case sym_enumeration_constant:
  case sym_INT: case sym_LONG: return "int";
  case sym_BOOL: return "bool";
  case sym_SCALAR: return "scalar";
  case sym_VECTOR: return "vector";
  case sym_TENSOR: return "tensor";
  case sym_COORD:  return "coord";
  case sym__COORD: return "_coord";
  case sym_VEC4:   return "vec4";
  case sym_IVEC:   return "ivec";
  }
  return "unknown_type";
}

trace
char * build_shader (External * externals, const ForeachData * loop,
		     const RegionParameters * region, const unsigned nwg[2])
{
  int Nl = region->level > 0 ? 1 << (region->level - 1) : N/Dimensions.x;
  char s[30];
  snprintf (s, 19, "%d", nconst > 0 ? nconst : 1);
  char a[20];
  snprintf (a, 19, "%g", nconst > 0 ? _constant[0] : 0);
  char * fs = str_append (NULL,
#if !_CUDA
                          "#version 430\n",
#endif
                          glsl_preproc,
			  "const int _nconst = ", s, ";\n"
			  "const real _constant[_nconst] = {", a);
  for (int i = 1; i < nconst; i++) {
    snprintf (a, 19, "%g", _constant[i]);
    fs = str_append (fs, ",", a);
  }
  fs = str_append (fs, "};\n"
#if _CUDA
                   "uniform real * _data;\n"
#else
		   "layout(std430, binding = 0)"
		   " restrict buffer _data_layout { real f[]; } _data"
#endif
                   );
  if (GPUContext.nssbo > 1) {
    snprintf (a, 19, "%d", GPUContext.nssbo);
    snprintf (s, 29, "%ld", GPUContext.max_ssbo_size/sizeof(real));
    fs = str_append (fs, "[", a, "];\n"
		     "#define _data_val(field,index) _data[_offset[(field)].i+(_offset[(field)].j+(index))/", s,
		     "].f[(_offset[(field)].j+(index))%", s, "]\n");
  }
  else
#if _CUDA
    fs = str_append (fs,
		     "#define _data_val(field,index) _data[(field)*field_size() + (index)]\n");
#else
    fs = str_append (fs, ";\n"
		     "#define _data_val(field,index) _data.f[(field)*field_size() + (index)]\n");
#endif
  
  /**
  Scalar field attributes */
      
  char * attributes = NULL;
  for (const External * g = externals; g; g = g->next)
    if (g->name[0] == '.') {
      attributes = str_append (attributes, "  ", type_string (g), " ", g->name + 1);
      for (int * d = g->data; d && *d > 0; d++) {
	char s[20]; snprintf (s, 19, "%d", *d);
	attributes = str_append (attributes, "[", s, "]");
      }
      attributes = str_append (attributes, ";\n");
    }
  
  if (attributes) {
    fs = str_append (fs, "struct _Attributes {\n", attributes, "};\n");
    sysfree (attributes);
    int nindex = 0;
    for (scalar s in baseblock)
      if (s.gpu.index)
	nindex++;
    assert (nindex > 0);
    char s[20]; snprintf (s, 19, "%d", nindex);
    fs = str_append (fs, "const _Attributes _attr[", s, "]={");
    nindex = 0;
    for (scalar s in baseblock)
      if (s.gpu.index) {
	fs = str_append (fs, nindex ? ",{" : "{"); nindex++;
	bool first = true;
	char * data = (char *) &_attribute[s.i];
	for (const External * g = externals; g; g = g->next)
	  if (g->name[0] == '.') {
	    if (!first) fs = str_append (fs, ",");
	    first = false;
	    if (g->type == sym_INT) {
	      char s[20]; snprintf (s, 19, "%d", *((int *)(data + g->nd)));
	      fs = str_append (fs, s);
	    }
	    else if (g->type == sym_BOOL)
	      fs = str_append (fs, *((bool *)(data + g->nd)) ? "true" : "false");
	    else if (g->type == sym_FLOAT) {
	      char s[20]; snprintf (s, 19, "%g", *((float *)(data + g->nd)));
	      fs = str_append (fs, s);
	    }
	    else if (g->type == sym_DOUBLE) {
	      char s[20]; snprintf (s, 19, "%g", *((double *)(data + g->nd)));
	      fs = str_append (fs, s);
	    }
	    else if (g->type == sym_IVEC) {
	      ivec * v = (ivec *)(data + g->nd);
	      char s[20]; snprintf (s, 19, "{%d,%d}", v->x, v->y);
	      fs = str_append (fs, s);
	    }
	    else if (g->type == sym_function_declaration) {
	      void * func = *((void **)(data + g->nd));
	      if (!func)
		fs = str_append (fs, "0");
	      else {
		External * ptr = _get_function ((long) func);
		char s[20]; snprintf (s, 19, "%d", ptr->nd);
		fs = str_append (fs, s);
	      }
	    }
	    else if (g->type == sym_SCALAR)
	      fs = write_scalar (fs, *((scalar *)(data + g->nd)));
	    else if (g->type == sym_VECTOR)
	      fs = write_vector (fs, *((vector *)(data + g->nd)));
	    else if (g->type == sym_TENSOR)
	      fs = write_tensor (fs, *((tensor *)(data + g->nd)));
	    else
	      fs = str_append (fs, "not implemented");
	  }
	fs = str_append (fs, "}");
      }
    fs = str_append (fs, "};\n");
  }

  /**
  Field offsets when using multiple SSBOs */

  if (GPUContext.nssbo > 1) {
    int imax = 1;
    for (scalar s in baseblock)
      if (s.gpu.index && s.i + s.block > imax)
	imax = s.i + s.block;
    for (External * g = externals; g; g = g->next)
      if (g->reduct) {
	scalar s = g->s;
	if (s.i + s.block > imax)
	  imax = s.i + s.block;
      }
    char string[100]; snprintf (string, 19, "%d", imax);
    fs = str_append (fs, "const struct { uint i, j; } _offset[", string, "]={");
    for (int s = 0; s < imax; s++) {
      fs = str_append (fs, s ? ",{" : "{");
      size_t offset = s*field_size(), size = GPUContext.max_ssbo_size/sizeof(real);
      size_t i = offset/size, j = offset%size;
      snprintf (string, 99, "%ld,%ld", i, j);
      fs = str_append (fs, string, "}");

      if (j + field_size() > 1L << 32) {
	fprintf (stderr, "%s:%d: error: resolution is too high for 32-bits indexing\n",
		 __FILE__, LINENO);
	exit (1);
      }
    }
    fs = str_append (fs, "};\n");
  }
  
  /**
  Non-local variables */

  for (External * g = externals; g; g = g->next) {
    if (g->name[0] == '.') {
      if (g->type == sym_function_declaration) {
	bool boundary = is_boundary_attribute (g);
	fs = str_append (fs, "#define _attr_", g->name + 1, "(s,args) (");
	foreach_function (f, f->used = false);
	char * expr = NULL;
	for (scalar s in baseblock)
	  if (s.gpu.index) {
	    char * data = (char *) &_attribute[s.i];
	    void * func = *((void **)(data + g->nd));
	    if (func) {
	      External * f = _get_function ((long) func);
	      if (!f->used && (!boundary || (s.stencil.io & s_output))) {
		f->used = true;
		char * args = is_void_function (f->data) ? " args,0" : " args";
		if (!expr)
		  expr = str_append (NULL, "(", f->name, args, ")");
		else {
		  char index[20];
		  snprintf (index, 19, "%d", f->nd);
		  char * s = str_append (NULL, "_attr(s,", g->name + 1, ")==", index,
					 "?(", f->name, args, "):", expr);
		  sysfree (expr);
		  expr = s;
		}
	      }
	    }
	  }
	if (expr) {
	  fs = str_append (fs, expr, ")\n");
	  sysfree (expr);
	}
	else
	  fs = str_append (fs, "0)\n");
      }
    }
    else if (g->type == sym_function_definition) {
      External * f = _get_function ((long) g->pointer);
#if _HIP_AMD
      fs = str_append (fs, "\n__device__ ", f->data, "\n");
#else
      fs = str_append (fs, "\n", f->data, "\n");
#endif
      char s[20]; snprintf (s, 19, "%d", f->nd);
      fs = str_append (fs, "const int _p", g->name, " = ", s, ";\n");
    }
    else if (g->type == sym_function_declaration) {
      External * f = _get_function ((long) g->pointer);
      char s[20]; snprintf (s, 19, "%d", f->nd);
      fs = str_append (fs, "const int ", g->name, " = ", s, ";\n",
		       "#define _f", g->name, "(args) (", f->name, " args)\n");
    }
    else if (g->type != sym_SCALAR &&
	     g->type != sym_VECTOR &&
	     g->type != sym_TENSOR) {
      if (g->type == sym_INT && !strcmp (g->name, "N")) {
	int level = region->level > 0 ? region->level - 1 : depth();
	char s[20], d[20], l[20];
	snprintf (s, 19, "%d", Nl);
	snprintf (d, 19, "%d", depth());
	snprintf (l, 19, "%d", level);
	fs = str_append (fs,
			 "const uint N = ", s, ", _depth = ", d, ";\n");
	if (GPUContext.nssbo <= 1) {
	  char size[30];
	  snprintf (size, 29, "%ld", (size_t) field_size());	  
	  fs = str_append (fs,
			   "const uint _field_size = ", size, ";\n");
	}
#ifdef shift_level // multigrid
	fs = str_append (fs,
			 "const uint _shift[_depth + 1] = {");
	for (int d = 0; d <= depth(); d++) {
	  snprintf (s, 19, "%ld", shift_level(d));
	  fs = str_append (fs,  d > 0 ? "," : "", s);
	}
	fs = str_append (fs, "};\n");
#endif // ifdef shift_level
        snprintf (s, 19, "%d", Dimensions.x);
        snprintf (d, 19, "%d", Dimensions.y);
	fs = str_append (fs,
                         "const ivec2 Dimensions = {", s, ",", d, "};\n"
			 "const uint NY = N*", d,
                         loop->face > 1 || loop->vertex ? "+1" : "", ";\n");
#if !_CUDA
	if (GPUContext.fragment_shader)
	  fs = str_append (fs, "in vec2 vsPoint;\n"
			   "Point point = {int((vsPoint.x*vsScale.x + vsOrigin.x)*N*Dimensions.x)"
			   " + GHOSTS,"
			   "int((vsPoint.y*vsScale.y + vsOrigin.y)*N*Dimensions.x) + GHOSTS,", l,
#if MULTIGRID
			   ",{int(N*Dimensions.x),int(N*Dimensions.y)}"
#else
			   ",N"
#endif
#if LAYERS
			   ",0"
#endif
			   "};\n"
			   "out vec4 FragColor;\n");
	else {
	  char nwgx[20], nwgy[20];
	  snprintf (nwgx, 19, "%d", nwg[0]);
	  snprintf (nwgy, 19, "%d", nwg[1]);
	  fs = str_append (fs, "layout (local_size_x = ", nwgx,
			   ", local_size_y = ", nwgy, ") in;\n");
	}
#endif // !_CUDA
      }
      else if (g->type == sym_INT && !strcmp (g->name, "nl")) {

	/**
	'int nl' gets special treatment. */
	
	char nl[20];
	snprintf (nl, 19, "%d", *((int *)g->pointer));
	fs = str_append (fs, "const int nl = ", nl, ";\n");
      }
      else if (g->type == sym_INT && !strcmp (g->name, "bc_period_x"))
	fs = str_append (fs, "const int bc_period_x = ", Period.x ?
			 "int(N*Dimensions.x)" : "-1", ";\n");
      else if (g->type == sym_INT && !strcmp (g->name, "bc_period_y"))
	fs = str_append (fs, "const int bc_period_y = ", Period.y ?
			 "int(N*Dimensions.y)" : "-1", ";\n");
      else if (GPUContext.fragment_shader && (region->n.x > 1 || region->n.y > 1) &&
	       g->type == sym_COORD && !strcmp (g->name, "p")) {

	/**
	'coord p' is assumed to be the parameter of a region. This is
	not flexible (the parameter must be called 'p') and should be improved. */
	
	fs = str_append (fs, "coord p = vec3((vsPoint*vsScale + vsOrigin)*L0 + vec2(X0, Y0),0);\n");
      }
      else if (IS_EXTERNAL_CONSTANT(g)) {
	// fixme: for the moment only 'const int' are considered, this could be generalised
	char value[20];
	assert (g->pointer);
	snprintf (value, 19, "%d", *((int *)g->pointer));
	fs = str_append (fs, "const ", type_string (g), " ", EXTERNAL_NAME (g), "=", value, ";\n");
      }
      else if (strcmp (g->name, "Dimensions")) {
	char * type = type_string (g);
	fs = str_append (fs, "uniform ", type, " ", EXTERNAL_NAME (g));
	for (int * d = g->data; d && *d > 0; d++) {
	  char s[20]; snprintf (s, 19, "%d", *d);
	  fs = str_append (fs, "[", s, "]");
	}
	fs = str_append (fs, ";\n");
      }
    }
    else { // scalar, vector and tensor fields
      int size = list_size (g);
      for (int j = 0; j < size; j++) {
	if (j == 0) {
	  fs = str_append (fs, "const ", type_string (g), " ", EXTERNAL_NAME (g));
	  if (g->nd == 0)
	    fs = str_append (fs, " = ");
	  else {
	    char s[20]; snprintf (s, 19, "%d", size);
	    fs = str_append (fs, "[", s, "] = {");
	  }
	}
	if (g->nd == 0 || j < size - 1) {
	  if (g->type == sym_SCALAR)
	    fs = write_scalar (fs, ((scalar *)g->pointer)[j]);
	  else if (g->type == sym_VECTOR)
	    fs = write_vector (fs, ((vector *)g->pointer)[j]);
	  else if (g->type == sym_TENSOR)
	    fs = write_tensor (fs, ((tensor *)g->pointer)[j]);
	  else
	    assert (false);
	}
	else { // last element of a list is always ignored (this is necessary for empty lists)
	  if (g->type == sym_SCALAR)
	    fs = str_append (fs, "{0,0}");
	  else if (g->type == sym_VECTOR)
	    fs = str_append (fs, "{{0,0},{0,0}}");
	  else if (g->type == sym_TENSOR)
	    fs = str_append (fs, "{{{0,0},{0,0}},{{0,0},{0,0}}}");
	  else
	    assert (false);
	}
	if (g->nd == 0)
	  fs = str_append (fs, ";\n");
	else if (j < size - 1)
	  fs = str_append (fs, ",");
	else
	  fs = str_append (fs, "};\n");
      }
    }
  }
    
  return fs;
}

trace
Shader * load_shader (const char * fs, uint32_t hash, const ForeachData * loop)
{
  assert (gpu_grid->shaders);
  khiter_t k = kh_get (INT, gpu_grid->shaders, hash);
  if (k != kh_end (gpu_grid->shaders)) {
    sysfree ((void *)fs);
    return kh_value (gpu_grid->shaders, k);
  }
#if PRINTSHADER
  {
    static int n = 1;
    fprintf (stderr, "=================== %s:%d: shader #%d ===================\n",
	     loop ? loop->fname : "reduction", loop ? loop->line : 0, n++);
    fputs (fs, stderr);
  }
#endif
  Shader * shader;
  if (loop && loop->func) {
    char func[strlen(loop->func) + 20];
    sprintf (func, "%s_%d", loop->func, loop->line);
    shader = load_normal_shader (fs, func, loop->fname, loop->line);
  }
  else
    shader = load_normal_shader (fs, __func__, __FILE__, LINENO);
  if (shader) {
    int ret;
    khiter_t k = kh_put (INT, gpu_grid->shaders, hash, &ret);
    assert (ret > 0);
    kh_value (gpu_grid->shaders, k) = shader;
  }
  sysfree ((void *)fs);
  return shader;
}

void gpu_free()
{
  if (!grid)
    return;
  free_boundaries();
  if (gpu_grid->shaders) {
    Shader * shader;
    int nshaders = 0;
    kh_foreach_value (gpu_grid->shaders, shader,
                      free_shader (shader);
                      nshaders++; );
#if PRINTNSHADERS  
    fprintf (stderr, "# %d shaders\n", nshaders);
#endif
    kh_destroy (INT, gpu_grid->shaders);
    gpu_grid->shaders = NULL;
    gpu_free_context (gpu_grid->data);
    gpu_grid->data = NULL;
  }
}

void gpu_free_grid (void)
{
  if (!grid)
    return;
  gpu_free();
  free_grid();
}

attribute {
  double (* boundary_left)   (Point, Point, scalar, bool *);
  double (* boundary_right)  (Point, Point, scalar, bool *);
  double (* boundary_top)    (Point, Point, scalar, bool *);
  double (* boundary_bottom) (Point, Point, scalar, bool *);
}

/**
The `stored` attibute tracks where the up-to-date field is stored:

*   0: on both the CPU and GPU (i.e. synchronized). 
*   1: on the CPU.
* - 1: on the GPU.
*/

attribute {
  struct {
    int stored, index;
  } gpu;
}

#define reset(...) reset_gpu (__VA_ARGS__)

trace
void reset_gpu (void * alist, double val)
{
  scalar * list = alist;
  for (scalar s in list)
    if (!is_constant(s)) {
      reset_scalar (s.i, s.block, field_size(), val);
      s.gpu.stored = -1;
    }
}

void gpu_init_grid (int n)
{
  if (grid && n == N)
    return;
  gpu_free_grid();
  init_grid (n);
  grid = realloc (grid, sizeof (GridGPU));
  /* GPUs drivers often generate floating-point exceptions... turn them off */
  disable_fpe (FE_DIVBYZERO|FE_INVALID);
  if (gpu_init_context (&gpu_grid->data)) {
    free_solver_func_add (gpu_free);
    free_solver_func_add (gpu_free_solver);
  }
  gpu_grid->shaders = kh_init (INT);
  realloc_ssbo (field_size());
}

// overload the default various functions

#define init_grid(n)  gpu_init_grid(n)
#undef  free_grid
#define free_grid()   gpu_free_grid()

static External * append_external (External * externals, External ** end, External * g)
{
  if (externals)
    (*end)->next = g;
  else
    externals = g;
  *end = g;
  (*end)->next = NULL;
  return externals;
}

static External * merge_external (External * externals, External ** end, External * g,
				  const ForeachData * loop)
{
  if (g->type == sym_function_declaration || g->type == sym_function_definition) {
    bool boundary = is_boundary_attribute (g);
    for (scalar s in baseblock)
      if (g->name[0] != '.' || (!boundary && s.gpu.index) ||
          (boundary && (s.stencil.io & s_output))) {
	void * ptr = g->name[0] != '.' ? g->pointer :
	  *((void **)(((char *) &_attribute[s.i]) + g->nd));
	if (ptr) {
	  External * p = _get_function ((long) ptr);
	  if (!p) {
	    fprintf (stderr, "%s:%d: GLSL: error: unregistered function pointer '%s'\n",
		     loop->fname, loop->line, g->name);
	    return NULL;
	  }
	  if (!p->used) {
	    p->used = true;
	    for (External * i = p->externals; i && i->name; i++)
	      externals = merge_external (externals, end, i, loop);
	    externals = append_external (externals, end, p);
	  }
	}
	if (g->name[0] != '.')
	    break;
      }
  }
  for (External * i = externals; i; i = i->next)
    if (!strcmp (g->name, i->name)) {

      /**
      Check whether a local *g* (resp. *i*) shadows a global *i* (resp
      *g*). */
      
      if (g->global == 0 && i->global == 1) g->global = 2;
      else if (g->global == 1 && i->global == 0) i->global = 2;
      else // already in the list
	return externals;
    }
  return append_external (externals, end, g);
}

static External * merge_externals (External * externals, const ForeachData * loop)
{
  External * merged = NULL, * end = NULL;
  static External ext[] = {
    { .name = "X0", .type = sym_DOUBLE, .pointer = &X0, .global = 1 },
    { .name = "Y0", .type = sym_DOUBLE, .pointer = &Y0, .global = 1 },
    { .name = "Z0", .type = sym_DOUBLE, .pointer = &Z0, .global = 1 },
    { .name = "L0", .type = sym_DOUBLE, .pointer = &L0, .global = 1 },
    { .name = "N",  .type = sym_INT,    .pointer = &N, .global = 1 },
#if LAYERS
    { .name = "nl",  .type = sym_INT, .pointer = &nl, .global = 1 },
    { .name = "_layer",  .type = sym_INT, .pointer = &_layer, .global = 1 },
    { .name = ".block", .type = sym_INT, .nd = attroffset (block) },
#endif
    { .name = NULL }
  };
  static External bc = {
    .name = "apply_bc", .type = sym_function_declaration, .pointer = (void *)(long)apply_bc
  };

  for (External * g = externals; g->name; g++) {
    g->used = false;
    if (g->global == 2) g->global = 0;
  }
  foreach_function (f, f->used = false);
  for (External * g = ext; g->name; g++) {
    g->used = false;
    merged = merge_external (merged, &end, g, loop);
  }
  if (loop->dirty) {
    bc.used = false;
    merged = merge_external (merged, &end, &bc, loop);
  }
  for (External * g = externals; g->name; g++)
    merged = merge_external (merged, &end, g, loop);
#if PRINTEXTERNAL  
  for (External * i = merged; i; i = i->next)
    fprintf (stderr, "external %s %d %p %d\n", i->name, i->type, i->pointer, i->global);
#endif
  if (loop->dirty)
    for (External * g = merged; g; g = g->next)
      if (g->global && !strcmp (g->name, "apply_bc_list"))
	g->pointer = loop->dirty;
  return merged;
}

static char * shader_append_func (char * s, const ForeachData * loop)
{
#if _CUDA
  char func[strlen(loop->func) + 20];
  sprintf (func, "%s_%d", loop->func, loop->line);
  return str_append (s,
                     "extern \"C\" __global__\n"
                     "void ", func, "() {\n");
#else // !_CUDA
  return str_append (s, "void main() {\n");
#endif
}

trace
static Shader * compile_shader (ForeachData * loop,
				uint32_t hash,
				const RegionParameters * region,
				External * externals,
				const char * kernel)
{  
  const char * error = strstr (kernel, "@error ");
  if (error) {
    for (const char * s = error + 7; *s != '\n' && *s != '\0'; s++)
      fputc (*s, stderr);
    fputc ('\n', stderr);
    loop->data = NULL;
    return NULL;
  }

  External * merged = merge_externals (externals, loop);
  if (!merged) {
    loop->data = NULL;
    return NULL;
  }
  
  int local = false;
  for (const External * g = merged; g; g = g->next) {
    if (g->global == 2)
      local = true;
    if (g->type != sym_SCALAR && g->type != sym_VECTOR && g->type != sym_TENSOR) {
      if (g->reduct && !strchr ("+mM", g->reduct)) {
	if (loop->first)
	  fprintf (stderr,
		   "%s:%d: GLSL: error: unknown reduction operation '%c'\n",
		   loop->fname, loop->line, g->reduct);
	return NULL;
      }
      if (g->type == sym_COORD || g->type == sym__COORD || g->type == sym_VEC4) {
	if (g->reduct) {
	  if (loop->first)
	    fprintf (stderr,
		     "%s:%d: GLSL: error: reductions not implemented for '%s' type\n",
		     loop->fname, loop->line, type_string (g));
	  return NULL;	
	}
      }
      else if (g->type != sym_FLOAT &&
	       g->type != sym_DOUBLE &&
	       g->type != sym_INT &&
	       g->type != sym_LONG &&
	       g->type != sym_BOOL &&
	       g->type != sym_enumeration_constant &&
	       g->type != sym_IVEC &&
	       g->type != sym_function_declaration &&
	       g->type != sym_function_definition) {
	if (loop->first)
	  fprintf (stderr, "%s:%d: GLSL: error: unknown type %d for '%s'\n",
		   loop->fname, loop->line, g->type, g->name);
	return NULL;
      }
    }
  }

  /**
  ## Number of compute shader work groups and groups */

  static const int NWG[2] = {16, 16};
  unsigned ng[2], nwg[2];
  int Nl = region->level > 0 ? 1 << (region->level - 1) : N/Dimensions.x;
  int * dims = &Dimensions.x;
  if (loop->face || loop->vertex)
    for (int i = 0; i < 2; i++) {
      if (Nl*dims[1-i] > NWG[i]) {
	nwg[i] = NWG[i] + 1;
	ng[i] = Nl*dims[1-i]/NWG[i];
      }
      else {
	nwg[i] = Nl*dims[1-i] + 1;
	ng[i] = 1;
      }
      assert (nwg[i]*ng[i] >= Nl*dims[1-i] + 1);
    }
  else
    for (int i = 0; i < 2; i++) {
      nwg[i] = Nl < NWG[i] ? Nl : NWG[i];
      ng[i] = Nl*dims[1-i]/nwg[i];
      assert (nwg[i]*ng[i] == Nl*dims[1-i]);
    }
 
  char * shader = build_shader (merged, loop, region, nwg);
  if (!shader)
    return NULL;

  /**
  ## main() */

  if (local) {
    shader = str_append (shader, "void _loop (");
    for (const External * g = merged; g; g = g->next)
      if (g->global == 2) {
	shader = str_append (shader, local++ == 1 ? "" : ", ", type_string (g), " ", g->name);
	if (g->nd) {
	  int size = list_size (g);
	  if (size > 0) {
	    char s[20]; snprintf (s, 19, "%d", size);
	    shader = str_append (shader, "[", s, "]");
	  }
	}
      }
    shader = str_append (shader, ") {\n");
  }
  else
    shader = shader_append_func (shader, loop);

  if (!GPUContext.fragment_shader) {
    char d[20];
    snprintf (d, 19, "%d", region->level > 0 ? region->level - 1 : depth());
    shader = str_append (shader,
#if _CUDA
                         "Point point = {csOrigin.x + int(blockIdx.y * blockDim.y + threadIdx.y) + GHOSTS,"
			 "csOrigin.y + int(blockIdx.x * blockDim.x + threadIdx.x) + GHOSTS,", d,
#if MULTIGRID
			 ",{(1<<",d,")*Dimensions.x,(1<<",d,")*Dimensions.y}"
#else
			 ",N"
#endif
#else // !_CUDA
                         "Point point = {csOrigin.x + int(gl_GlobalInvocationID.y) + GHOSTS,"
			 "csOrigin.y + int(gl_GlobalInvocationID.x) + GHOSTS,", d,
#if MULTIGRID
			 ",{(1<<",d,")*Dimensions.x,(1<<",d,")*Dimensions.y}"
#else
			 ",N"
#endif
#endif // !_CUDA
#if LAYERS
			 ",0};\n"
#else
			 "};\n"
#endif
			 );
  }
  shader = str_append (shader,
		       "if (point.i < N*Dimensions.x + 2*GHOSTS && "
		       "point.j < N*Dimensions.y + 2*GHOSTS) {\n");
  if (loop->vertex)
    shader = str_append (shader, "  int ig = -1, jg = -1;\n");
#if 1
  for (const External * g = merged; g; g = g->next)
    if (g->reduct) {
      shader = str_append (shader, type_string (g), " ",
                           g->global == 2 ? "_loc_" : "", g->name, " = ",
                           EXTERNAL_NAME (g), ";\n");
      shader = str_append (shader, "const scalar ", g->name, "_out_ = ");
      shader = write_scalar (shader, g->s);
      shader = str_append (shader, ";\n");
    }
#endif
  
  shader = str_append (shader, kernel);
  shader = str_append (shader, "\nif (point.j - GHOSTS < NY) {");
  for (const External * g = merged; g; g = g->next)
    if (g->reduct) {
      shader = str_append (shader, "\n  val_red_(", g->name, "_out_) = ", g->name, ";");
      scalar s = g->s;
      s.gpu.stored = -1;
    }
  shader = str_append (shader, "\n}",
		       loop->dirty ? "apply_bc(point);" : "",
		       "}}\n");

  if (local) {
    shader = shader_append_func (shader, loop);
    shader = str_append (shader, "_loop(");
    local = 1;
    for (const External * g = merged; g; g = g->next)
      if (g->global == 2)
	shader = str_append (shader, local++ == 1 ? "" : ",", EXTERNAL_NAME (g));
    shader = str_append (shader, ");}\n");
  }
    
  Shader * s = load_shader (shader, hash, loop);
  loop->data = s;
  if (!s)
    return NULL;

  finalize_shader (s, externals, merged, ng, nwg);
  
  return s;
}

static
void free_reduction_fields (const External * externals)
{
  for (const External * g = externals; g; g = g->next)
    if (g->reduct) {
      scalar s = g->s;
      delete ({s});
    }
}

trace
static void gpu_cpu_sync (scalar * list, SyncMode mode, const char * fname, int line)
{
#if PRINTCOPYGPU
  bool copy = false;
#endif
  for (scalar s in list)
    if (((s.stencil.io & s_input) || (s.stencil.io & s_output)) &&
        ((mode == GPU_READ && s.gpu.stored < 0) ||
         (mode == GPU_WRITE && s.gpu.stored > 0))) {
      if (s.gpu.stored > 0 && !(s.stencil.bc & s_centered))
        boundary ({s});
#if PRINTCOPYGPU
      if (!copy) {
	fprintf (stderr, "%s:%d: %s ", fname, line,
		 mode == GPU_READ ? "importing" : "exporting");
	copy = true;
	gpu_cpu_sync_scalar (s.i, s.block, grid_data(), mode);
        fprintf (stderr, "{%s", s.name);
      }
      else {
	gpu_cpu_sync_scalar (s.i, s.block, grid_data(), mode);
        fprintf (stderr, ",%s", s.name);
      }
#else
      gpu_cpu_sync_scalar (s.i, s.block, grid_data(), field_size(), mode);
#endif
      s.gpu.stored = 0;
    }
#if PRINTCOPYGPU
  if (copy)
    fprintf (stderr, "} %s GPU\n", mode == GPU_READ ? "from" : "to");
#endif
}

trace
static Shader * setup_shader (ForeachData * loop, const RegionParameters * region,
			      External * externals,
			      const char * kernel)
{
  /**
  We will directly apply boundary conditions to fields marked 'dirty'
  by automatic stencils. */
  
  apply_bc_list = loop->dirty;
  for (scalar s in loop->dirty) {
    
    /**
    We also apply boundary stencils so that input/output are also set
    properly for boundary conditions which may use external fields. */
  
    for (int b = 0; b < nboundary; b++)
      if (s.boundary_stencil[b])
        s.boundary_stencil[b] ((Point){0}, (Point){0}, s, NULL);

    /**
    We make sure all fields marked dirty are also outputs. */

    if (!(s.stencil.io & s_output))
      s.stencil.io |= s_output;
  }

  for (scalar s in baseblock)
    s.gpu.index = 0;
  int index = 1;
  for (scalar s in baseblock)
    if (((s.stencil.io & s_input) || (s.stencil.io & s_output)) && !s.gpu.index) {
#if PRINTBOUNDARY
      fprintf (stderr, "%s:%d: %s:%s%s index: %d\n",
               loop->fname, loop->line, 
               s.name, (s.stencil.io & s_input) ? " input" : "",
               (s.stencil.io & s_output) ? " output" : "", index);
#endif
      if (s.v.x.i == -1) // scalar
	s.gpu.index = index++;
      else { // vector
	vector v = s.v;
	for (scalar c in {v})
	  if (!c.gpu.index)
	    c.gpu.index = index++;
      }
    }
  for (scalar s in loop->dirty) {
#if PRINTBOUNDARY
    fprintf (stderr, "%s:%d: dirty: %s:%s%s index: %d\n",
             loop->fname, loop->line, s.name,
             (s.stencil.io & s_input) ? " input" : "",
             (s.stencil.io & s_output) ? " output" : "",
             s.gpu.index);
#endif
    s.boundary_left   = s.boundary[left];
    s.boundary_right  = s.boundary[right];
    s.boundary_top    = s.boundary[top];
    s.boundary_bottom = s.boundary[bottom];
  }
  
  /**
  ## Allocate reduction fields */

  for (External * g = externals; g && g->name; g++)
    if (g->reduct) {
      scalar s = new scalar;
      s.stencil.io |= s_output;
      g->s = s;
#if PRINTREDUCT
      fprintf (stderr, "%s:%d: new reduction field %d for %s\n",
	       loop->fname, loop->line, s.i, g->name);
#endif
    }

  /**
  ## Reuse or compile a new shader */
  
  Shader * shader;
  uint32_t hash = hash_shader (externals, loop, region, kernel);
  assert (gpu_grid->shaders);
  khiter_t k = kh_get (INT, gpu_grid->shaders, hash);
  if (k != kh_end (gpu_grid->shaders))
    shader = kh_value (gpu_grid->shaders, k);
  else {
    shader = compile_shader (loop, hash, region, externals, kernel);  
    if (!shader) {
      free_reduction_fields (externals);
      return NULL;
    }
  }

  gpu_cpu_sync (baseblock, GPU_WRITE, loop->fname, loop->line);

  /**
  ## Apply boundary conditions
  
  This can be required if boundary conditions have been modified
  between loops. */

  scalar * listc = NULL;
  for (scalar s in loop->listc)
    if (!(s.stencil.bc & s_centered))
      listc = list_prepend (listc, s);
  scalar * listf_x = NULL, * listf_y = NULL;
  foreach_dimension()
    for (scalar s in loop->listf.x)
      if (!(s.stencil.bc & s_face))
	listf_x = list_prepend (listf_x, s);
  if (listc || listf_x || listf_y) {
#if PRINTBC
    fprintf (stderr, "%s:%d: applying BCs for", loop->fname, loop->line);
    for (scalar s in listc)
      fprintf (stderr, " %s", s.name);
    foreach_dimension()
      for (scalar s in listf_x)
	fprintf (stderr, " %s", s.name);
    fputc ('\n', stderr);
#endif
    int nvar = datasize/sizeof(real);
    _Attributes backup[nvar];
    memcpy (backup, _attribute, nvar*sizeof(_Attributes));
    foreach (gpu) {
      for (scalar s in listc)
	s[] = s[]; // does nothing but ensures that BCs are applied at the end of the loop
      foreach_dimension()
	for (scalar s in listf_x)
	  s[] = s[];
    }
    memcpy (_attribute, backup, nvar*sizeof(_Attributes));
    for (scalar s in listc) {
      s.stencil.bc |= s_centered;
      s.gpu.stored = -1;
    }
    free (listc); 
    foreach_dimension() {
      for (scalar s in listf_x) {
        s.stencil.bc |= s_face;
	s.gpu.stored = -1;
      }
      free (listf_x);
    }
    apply_bc_list = loop->dirty;
  }
  
  post_setup_shader (shader, externals);
  
  return shader;
}

static bool doloop_on_gpu (ForeachData * loop, const RegionParameters * region,
			   External * externals,
			   const char * kernel)
{
  Shader * shader = setup_shader (loop, region, externals, kernel);
  if (!shader)
    return false;
  
  /**
  ## Render */

  int Nl = run_shader (shader, region);

  /**
  ## Perform reductions and cleanup */

  bool nreductions = false;
  for (const External * g = externals; g && g->name; g++)
    if (g->reduct) {
      nreductions = true;
      break;
    }
  if (nreductions)
    tracing ("gpu_reduction", loop->fname, loop->line);
  for (const External * g = externals; g && g->name; g++)
    if (g->reduct) {
      scalar s = g->s;
      double result = gpu_reduction (field_offset(s, region->level), g->reduct, region,
                                     gpu_grid->data,
				     loop->face == 1 ?  (Nl*Dimensions.x + 1)*Nl*Dimensions.y :
				     loop->face == 2 ?   Nl*Dimensions.x*(Nl*Dimensions.y + 1) :
				     loop->face == 3 || loop->vertex ?
				     (Nl*Dimensions.x + 1)*(Nl*Dimensions.y + 1) :
				     sq(Nl)*Dimensions.x*Dimensions.y);
#if PRINTREDUCT
      fprintf (stderr, "%s:%d: %s %c %g\n",
	       loop->fname, loop->line, g->name, g->reduct, result);
#endif
      if (g->type == sym_DOUBLE) *((double *)g->pointer) = result;
      else if (g->type == sym_FLOAT) *((float *)g->pointer) = result;
      else if (g->type == sym_INT) *((int *)g->pointer) = result;
      else
	assert (false);
      delete ({s});
    }
  if (nreductions)
    end_tracing ("gpu_reduction", loop->fname, loop->line);

  return true;
}

bool gpu_end_stencil (ForeachData * loop,
		      const RegionParameters * region,
		      External * externals,
		      const char * kernel)
{
  bool on_gpu = ((loop->parallel == 1 && !on_cpu) || loop->parallel == 3) &&
    (loop->first || loop->data);
  if (on_gpu) {
    on_gpu = doloop_on_gpu (loop, region, externals, kernel);
    if (!on_gpu) {
      fprintf (stderr, "%s:%d: %s: foreach() done on CPU (see GLSL errors above)\n",
	       loop->fname, loop->line, loop->parallel == 3 ? "error" : "warning");
      if (loop->parallel == 3) // must run on GPU but cannot run
	exit (1);
      loop->data = NULL;
    }
  }
  if (on_gpu) {
    // do not apply BCs on CPU
    free (loop->listc), loop->listc = NULL;
    foreach_dimension()
      free (loop->listf.x), loop->listf.x = NULL;
    for (scalar s in loop->dirty)
      s.stencil.bc = s_centered|s_face;
    free (loop->dirty), loop->dirty = NULL;
  }
  else {
    gpu_cpu_sync (baseblock, GPU_READ, loop->fname, loop->line);
    boundary_stencil (loop);
  }
  
  for (scalar s in baseblock)
    if (s.stencil.io & s_output)
      s.gpu.stored = on_gpu ? -1 : 1;

  return on_gpu && loop->parallel != 3;
}

/**
## Useful (and not so useful) links

* [Core Language (GLSL)](https://www.khronos.org/opengl/wiki/Core_Language_(GLSL))
* [Image Load Store](https://www.khronos.org/opengl/wiki/Image_Load_Store)
* [Persistent Mapped Buffers in OpenGL](https://www.cppstories.com/2015/01/persistent-mapped-buffers-in-opengl/)
* [Best Practices for Working with Texture Data (OpenGL Programming Guide for Mac)](https://developer.apple.com/library/archive/documentation/GraphicsImaging/Conceptual/OpenGL-MacProgGuide/opengl_texturedata/opengl_texturedata.html)
* [https://github.com/Erkaman/fluid_sim]()
* [https://stackoverflow.com/questions/67529148/how-can-i-optimize-a-compute-shader-thats-heavy-on-imageload-calls]()
* [https://www.slideshare.net/Mark_Kilgard/siggraph-2012-nvidia-opengl-for-2012]()
* [https://stackoverflow.com/questions/7954927/passing-a-list-of-values-to-fragment-shader]()
* [https://computergraphics.stackexchange.com/questions/5416/why-use-image-load-stores-instead-of-fbos]()
* [https://computergraphics.stackexchange.com/questions/9956/performance-of-compute-shaders-vs-fragment-shaders-for-deferred-rendering]()
* [https://computergraphics.stackexchange.com/questions/54/when-is-a-compute-shader-more-efficient-than-a-pixel-shader-for-image-filtering]()
* [http://diaryofagraphicsprogrammer.blogspot.com/2014/03/compute-shader-optimizations-for-amd.html]()
* [DirectCompute: Optimizations and Best Practices, Eric Young, NVIDIA Corporation, 2010](http://on-demand.gputechconf.com/gtc/2010/presentations/S12312-DirectCompute-Pre-Conference-Tutorial.pdf)
* [Compute shaders in graphics: Gaussian blur](https://lisyarus.github.io/blog/graphics/2022/04/21/compute-blur.html)
* [Arm Mali GPUs Best Practices Developer Guide](https://armkeil.blob.core.windows.net/developer/Arm%20Developer%20Community/PDF/Arm%20Mali%20GPU%20Best%20Practices.pdf)
* [Rendering to a 3D texture](https://community.khronos.org/t/rendering-to-a-3d-texture/75285/2)
* [GLFW Shared Contexts](https://www.glfw.org/docs/3.3/context_guide.html#context_sharing)
* [Slides on using Array or Bindless Textures](https://www.slideshare.net/CassEveritt/beyond-porting)
* [Optimizing Compute Shaders for L2 Locality using Thread-Group ID Swizzling](https://developer.nvidia.com/blog/optimizing-compute-shaders-for-l2-locality-using-thread-group-id-swizzling/)
* [Optimizing GPU occupancy and resource usage with large thread groups](https://gpuopen.com/learn/optimizing-gpu-occupancy-resource-usage-large-thread-groups/)
* [https://www.reddit.com/r/CUDA/comments/lkhcbv/is_there_a_list_of_gpus_ranked_by_fp64/]()
* [https://www.techpowerup.com/gpu-specs/]()
*/
