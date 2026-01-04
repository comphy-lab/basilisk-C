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

## Bugs

* Trying to use more than one shader storage buffer (SSBO) with Intel
  drivers does not seem work. This limits the amount of available
  video memory to a single SSBO with a maximum size which is usually
  2GB (or less).

## See also

* [Test cases on GPUs](/src/test/README#gpu-tests)
* [GPU benchmarks](Benchmarks.md)
* [Computation kernels](/src/ast/kernels.c)

# Implementation */

bool on_cpu = false;

#include <khash.h>

typedef struct {
  int location, type, nd, local;
  void * pointer;
} MyUniform;

typedef struct {
  GLuint id, ng[2];
  MyUniform * uniforms;
} Shader;

KHASH_MAP_INIT_INT(INT, Shader *)
  
typedef struct {
  GRIDPARENT parent;
  GLuint reduct[2];
  khash_t(INT) * shaders;
} GridGPU;

#define gpu_grid ((GridGPU *)grid)

static char * str_append_array (char * dst, const char * list[])
{
  int empty = (dst == NULL);
  int len = empty ? 0 : strlen (dst);
  for (const char ** s = list; *s != NULL; s++)
    len += strlen (*s);
  dst = (char *) sysrealloc (dst, len + 1);
  if (empty) dst[0] = '\0';
  for (const char ** s = list; *s != NULL; s++)
    strcat (dst, *s);
  return dst;
}

#define str_append(dst, ...) str_append_array (dst, (const char *[]){__VA_ARGS__, NULL})
#define xstr(a) str(a)
#define str(a) #a

static char glsl_preproc[] =
  "// #line " xstr(LINENO) " " __FILE__ "\n"
  "#define dimensional(x)\n"
  "#define qassert(file, line, cond)\n"
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
  "#define _coord dvec2\n"
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
  "#define forin(type,s,list) for (int _i = 0; _i < list.length() - 1; _i++) { type s = list[_i];\n"
  "#define endforin() }\n"
#if LAYERS
  "#define _index(a,m) ((a).i + (point.l + _layer + (m) < _attr(a,block) ? point.l + _layer + (m) : 0))\n"
#else
  "#define _index(a,m) ((a).i)\n"
#endif
  "#define forin2(a,b,c,d) for (int _i = 0; _i < c.length() - 1; _i++)"
  "  { a = c[_i]; b = d[_i];\n"
  "#define endforin2() }\n"
  "#define forin3(a,b,e,c,d,f) for (int _i = 0; _i < c.length() - 1; _i++)"
  "  { a = c[_i]; b = d[_i]; e = f[_i];\n"
  "#define endforin3() }\n"
  "#define NOT_UNUSED(x)\n"
  "#define pi 3.14159265359\n"
  "#define nodata (1e30)\n"
  "#define fabs(x) abs(x)\n"
  "#define neighborp(_i,_j,_k) Point(point.i+_i,point.j+_j,point.level,point.n"
#if LAYERS
  ",point.l"
#endif
  ")\n"
  "const real z = 0.;\n"
  "const int ig = 0, jg = 0;\n"
  "layout (location = 0) uniform ivec2 csOrigin = ivec2(0,0);\n"
  "layout (location = 1) uniform vec2 vsOrigin = vec2(0.,0.);\n"
  "layout (location = 2) uniform vec2 vsScale = vec2(1.,1.);\n"
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
      if (g->name[0] != '.' || (!boundary && s.gpu.index) || (boundary && s.output)) {
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

#define EXTERNAL_NAME(g) (g)->global == 2 ? "_loc_" : "", (g)->name, (g)->reduct ? "_in_" : ""

trace
char * build_shader (External * externals, const ForeachData * loop,
		     const RegionParameters * region, const GLuint nwg[2])
{
  int Nl = region->level > 0 ? 1 << (region->level - 1) : N/Dimensions.x;
  char s[30];
  snprintf (s, 19, "%d", nconst > 0 ? nconst : 1);
  char a[20];
  snprintf (a, 19, "%g", nconst > 0 ? _constant[0] : 0);
  char * fs = str_append (NULL, "#version 430\n", glsl_preproc,
			  "const int _nconst = ", s, ";\n"
			  "const real _constant[_nconst] = {", a);
  for (int i = 1; i < nconst; i++) {
    snprintf (a, 19, "%g", _constant[i]);
    fs = str_append (fs, ",", a);
  }
  fs = str_append (fs, "};\n"
		   "layout(std430, binding = 0)"
		   " restrict buffer _data_layout { real f[]; } _data");
  if (GPUContext.nssbo > 1) {
    snprintf (a, 19, "%d", GPUContext.nssbo);
    snprintf (s, 29, "%ld", GPUContext.max_ssbo_size/sizeof(real));
    fs = str_append (fs, "[", a, "];\n"
		     "#define _data_val(field,index) _data[_offset[(field)].i+(_offset[(field)].j+(index))/", s,
		     "].f[(_offset[(field)].j+(index))%", s, "]\n");
  }
  else
    fs = str_append (fs, ";\n"
		     "#define _data_val(field,index) _data.f[(field)*field_size() + (index)]\n");
  
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
	      if (!f->used && (!boundary || s.output)) {
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
      fs = str_append (fs, "\n", f->data, "\n");
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
	if (GPUContext.nssbo == 1) {
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
	snprintf (s, 19, "{%d,%d}", Dimensions.x, Dimensions.y);
	fs = str_append (fs,
			 "const ivec Dimensions = ", s, ";\n"
			 "const uint NY = ", loop->face > 1 || loop->vertex ?
			 "N*Dimensions.y + 1" : "N*Dimensions.y", ";\n");
	if (GPUContext.fragment_shader)
	  fs = str_append (fs, "in vec2 vsPoint;\n"
			   "Point point = {int((vsPoint.x*vsScale.x + vsOrigin.x)*N*Dimensions.x)"
			   " + GHOSTS,"
			   "int((vsPoint.y*vsScale.y + vsOrigin.y)*N*Dimensions.x) + GHOSTS,", l,
#if MULTIGRID
			   ",ivec2(N*Dimensions.x,N*Dimensions.y)"
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
      else if (strcmp (g->name, "Dimensions")) {
	char * type = type_string (g);
	fs = str_append (fs, "uniform ", type, " ", EXTERNAL_NAME (g));
	for (int * d = g->data; d && *d > 0; d++) {
	  char s[20]; snprintf (s, 19, "%d", *d);
	  fs = str_append (fs, "[", s, "]");
	}
	fs = str_append (fs, ";\n");
	if (g->reduct) {
	  fs = str_append (fs, type, " ", g->global == 2 ? "_loc_" : "", g->name, " = ",
			   EXTERNAL_NAME (g), ";\n");
	  fs = str_append (fs, "const scalar ", g->name, "_out_ = ");
	  fs = write_scalar (fs, g->s);
	  fs = str_append (fs, ";\n");
	}
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
  GLuint id;
  if (!GPUContext.fragment_shader)
    id = loadNormalShader (NULL, fs);
  else {
    const char quad[] =
      "#version 430\n"
      "layout(location = 0) in vec3 vsPos;"
      "out vec2 vsPoint;"
      "void main() {"
      "  vsPoint = vsPos.xy;"
      "  gl_Position =  vec4(2.*vsPos.xy - vec2(1.), 0., 1.);"
      "}";
    id = loadNormalShader (quad, fs);
  }
  Shader * shader = NULL;
  if (id) {
    shader = calloc (1, sizeof (Shader));
    shader->id = id;
    int ret;
    khiter_t k = kh_put (INT, gpu_grid->shaders, hash, &ret);
    assert (ret > 0);
    kh_value (gpu_grid->shaders, k) = shader;
  }
  sysfree ((void *)fs);
  return shader;
}

void gpu_limits (FILE * fp)
{
  GLString * i = gpu_limits_list;
  while (i->s) {
    GLint val;
    GL_C (glGetIntegerv (i->index, &val));  
    fprintf (fp, "%s: %d\n", i->s, val);
    i++;
  }
}

void gpu_free()
{
  if (!grid)
    return;
  free_boundaries();
  Shader * shader;
  int nshaders = 0;
  kh_foreach_value (gpu_grid->shaders, shader,
		    free (shader->uniforms);
		    free (shader);
		    nshaders++; );
#if PRINTNSHADERS  
  fprintf (stderr, "# %d shaders\n", nshaders);
#endif
  kh_destroy (INT, gpu_grid->shaders);
  gpu_grid->shaders = NULL;
  if (gpu_grid->reduct[0]) {
    GL_C (glDeleteBuffers (2, gpu_grid->reduct));
    for (int i = 0; i < 2; i++)
      gpu_grid->reduct[i] = 0;
  }
  if (GPUContext.nssbo) {
    GL_C (glDeleteBuffers (GPUContext.nssbo, GPUContext.ssbo));
    free (GPUContext.ssbo);
    GPUContext.ssbo = NULL;
    GPUContext.nssbo = 0;
  }
  GPUContext.current_size = 0;
}

void gpu_init()
{
  /* GPUs drivers often generate floating-point exceptions... turn them off */
  disable_fpe (FE_DIVBYZERO|FE_INVALID);
  if (!GPUContext.window) {
    if (!glfwInit ())
      exit (1);

    glfwWindowHint (GLFW_VISIBLE, GL_FALSE);
    glfwWindowHint (GLFW_RESIZABLE, GL_FALSE);
    glfwWindowHint (GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint (GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint (GLFW_SAMPLES, 0);
    glfwWindowHint (GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint (GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#if DEBUG_OPENGL    
    glfwWindowHint (GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE);
#endif
    
    GPUContext.window = glfwCreateWindow (1, 1, "GPU grid", NULL, NULL);
    if (!GPUContext.window) {
      glfwTerminate();
      fprintf (stderr, "GLFW: error: could not create window!\n");
      exit (1);
    }
    glfwMakeContextCurrent (GPUContext.window);

    // load GLAD.
    assert (gladLoadGLLoader ((GLADloadproc)glfwGetProcAddress));
    assert (glBindImageTexture);

#if DEBUG_OPENGL    
    GLint flags;
    GL_C (glGetIntegerv (GL_CONTEXT_FLAGS, &flags));
    if (flags & GL_CONTEXT_FLAG_DEBUG_BIT) {
      GL_C (glEnable(GL_DEBUG_OUTPUT));
      GL_C (glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS));
      GL_C (glDebugMessageCallback (GLDebugMessageCallback, NULL));
      GL_C (glDebugMessageControl (GL_DONT_CARE, GL_DONT_CARE, GL_DONT_CARE, 0, NULL, GL_TRUE));
    }
#endif // DEBUG_OPENGL

    free_solver_func_add (gpu_free);
    free_solver_func_add (gpu_free_solver);
  }
  gpu_grid->shaders = kh_init (INT);
  for (int i = 0; i < 2; i++)
    gpu_grid->reduct[i] = 0;
    
  realloc_ssbo();
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

trace
static void gpu_cpu_sync_scalar (scalar s, char * sep, GLenum mode)
{
  assert ((mode == GL_MAP_READ_BIT && s.gpu.stored < 0) ||
	  (mode == GL_MAP_WRITE_BIT && s.gpu.stored > 0));
  if (s.gpu.stored > 0 && s.dirty)
    boundary ({s});
  GL_C (glMemoryBarrier (GL_BUFFER_UPDATE_BARRIER_BIT));
  size_t size = (size_t)field_size()*sizeof(real), offset = s.i*size, totalsize = s.block*size;
  char * cd = grid_data() + offset;
  int index = offset/GPUContext.max_ssbo_size;
  offset -= index*GPUContext.max_ssbo_size;
  while (totalsize) {
    GL_C (glBindBuffer (GL_SHADER_STORAGE_BUFFER, GPUContext.ssbo[index]));
    size_t size = min (totalsize, GPUContext.max_ssbo_size - offset);

    //    fprintf (stderr, "map %d %ld %ld\n", index, offset, size);
    
    char * gd = glMapBufferRange (GL_SHADER_STORAGE_BUFFER, offset, size, mode);
    assert (gd);
    if (mode == GL_MAP_READ_BIT)
      memcpy (cd, gd, size);
    else if (mode == GL_MAP_WRITE_BIT)
      memcpy (gd, cd, size);
    else
      assert (false);
    assert (glUnmapBuffer (GL_SHADER_STORAGE_BUFFER));
    cd += size, totalsize -= size, offset = 0, index++;
  }
  GL_C (glBindBuffer (GL_SHADER_STORAGE_BUFFER, 0));
  if (sep)
    fprintf (stderr, "%s%s", sep, s.name);
  s.gpu.stored = 0;
}

static void gpu_cpu_sync (scalar * list, GLenum mode, const char * fname, int line)
{
#if PRINTCOPYGPU
  bool copy = false;
#endif
  for (scalar s in list)
    if ((s.input || s.output) && ((mode == GL_MAP_READ_BIT && s.gpu.stored < 0) ||
				  (mode == GL_MAP_WRITE_BIT && s.gpu.stored > 0))) {
#if PRINTCOPYGPU
      if (!copy) {
	fprintf (stderr, "%s:%d: %s ", fname, line,
		 mode == GL_MAP_READ_BIT ? "importing" : "exporting");
	copy = true;
	gpu_cpu_sync_scalar (s, "{", mode);
      }
      else
	gpu_cpu_sync_scalar (s, ",", mode);
#else
      gpu_cpu_sync_scalar (s, NULL, mode);
#endif
    }
#if PRINTCOPYGPU
  if (copy)
    fprintf (stderr, "} %s GPU\n", mode == GL_MAP_READ_BIT ? "from" : "to");
#endif
}

trace
void reset_gpu (void * alist, double val)
{
  size_t size = (size_t)field_size()*sizeof(real);
  scalar * list = alist;
  for (scalar s in list)
    if (!is_constant(s)) {
      size_t offset = s.i*size, totalsize = max(s.block, 1)*size;
      int index = offset/GPUContext.max_ssbo_size;
      offset -= index*GPUContext.max_ssbo_size;
      while (totalsize) {
	GL_C (glBindBuffer (GL_SHADER_STORAGE_BUFFER, GPUContext.ssbo[index]));
	size_t size = min (totalsize, GPUContext.max_ssbo_size - offset);
#if SINGLE_PRECISION
	float fval = val;
	GL_C (glClearBufferSubData (GL_SHADER_STORAGE_BUFFER, GL_R32F,
				    offset, size,
				    GL_RED, GL_FLOAT, &fval));
#else
	GL_C (glClearBufferSubData (GL_SHADER_STORAGE_BUFFER, GL_RG32UI,
				    offset, size,
				    GL_RG_INTEGER, GL_UNSIGNED_INT, &val));      
#endif
	totalsize -= size, offset = 0, index++;
      }
      s.gpu.stored = -1;
    }
  GL_C (glBindBuffer (GL_SHADER_STORAGE_BUFFER, 0));
}

#define reset(...) reset_gpu (__VA_ARGS__)

void gpu_init_grid (int n)
{
  if (grid && n == N)
    return;
  gpu_free_grid();
  init_grid (n);
  grid = realloc (grid, sizeof (GridGPU));
  gpu_init();
}

// overload the default various functions

#define init_grid(n)  gpu_init_grid(n)
#undef  free_grid
#define free_grid()   gpu_free_grid()

#include "reduction.h"

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
      if (g->name[0] != '.' || (!boundary && s.gpu.index) || (boundary && s.output)) {
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
  GLuint ng[2], nwg[2];
  int Nl = region->level > 0 ? 1 << (region->level - 1) : N/Dimensions.x;
  int * dims = &Dimensions.x;
  if (loop->face || loop->vertex) {
    for (int i = 0; i < 2; i++)
      if (Nl*dims[1-i] > NWG[i]) {
	nwg[i] = NWG[i] + 1;
	ng[i] = Nl*dims[1-i]/NWG[i];
	assert (nwg[i]*ng[i] >= Nl*dims[1-i] + 1);
      }
      else {
	nwg[i] = Nl*dims[1-i] + 1;
	ng[i] = 1;
      }
  }
  else
    for (int i = 0; i < 2; i++) {
      nwg[i] = Nl*dims[1-i] > NWG[i] ? NWG[i] : Nl*dims[1-i];
      ng[i] = Nl*dims[1-i]/nwg[i];
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
    shader = str_append (shader, "void main() {\n");
  
  if (!GPUContext.fragment_shader) {
    char d[20];
    snprintf (d, 19, "%d", region->level > 0 ? region->level - 1 : depth());
    shader = str_append (shader, "Point point = {csOrigin.x + int(gl_GlobalInvocationID.y) + GHOSTS,"
			 "csOrigin.y + int(gl_GlobalInvocationID.x) + GHOSTS,", d,
#if MULTIGRID
			 ",ivec2((1<<",d,")*Dimensions.x,(1<<",d,")*Dimensions.y)"
#else
			 ",N"
#endif
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
    shader = str_append (shader, "void main(){_loop(");
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
  s->ng[0] = ng[0], s->ng[1] = ng[1];
  
  /**
  ## Make list of uniforms */

  for (External * g = merged; g; g = g->next)
    g->used = 0;
  int index = 1;
  for (External * g = externals; g && g->name; g++)
    g->used = index++;
  int nuniforms = 0;
  for (const External * g = merged; g; g = g->next) {
    if (g->name[0] == '.') continue;
    if (g->type == sym_function_declaration || g->type == sym_function_definition) continue;
    if (g->type == sym_INT && (!strcmp (g->name, "N") ||
			       !strcmp (g->name, "nl") ||
			       !strcmp (g->name, "bc_period_x") ||
			       !strcmp (g->name, "bc_period_y")))
      continue;
    if (g->type == sym_INT ||
	g->type == sym_LONG ||
	g->type == sym_FLOAT ||
	g->type == sym_DOUBLE ||
	g->type == sym__COORD ||
	g->type == sym_COORD ||
	g->type == sym_BOOL ||
	g->type == sym_VEC4) {
      char * name = str_append (NULL, EXTERNAL_NAME (g));
      int location = glGetUniformLocation (s->id, name);
      sysfree (name);
      if (location >= 0) {
	// fprintf (stderr, "%s:%d: %s\n", loop->fname, loop->line, name);
	// not an array or just a one-dimensional array
	assert (!g->nd);
	assert (!g->data || ((int *)g->data)[1] == 0);
	int nd = g->data ? ((int *)g->data)[0] : 1;
	s->uniforms = realloc (s->uniforms, (nuniforms + 2)*sizeof(MyUniform));
	s->uniforms[nuniforms] = (MyUniform){
	  .location = location, .type = g->type, .nd = nd,
	  .local = g->global == 1 ? -1 : g->used - 1,
	  .pointer = g->global == 1 ? g->pointer : NULL };
	s->uniforms[nuniforms + 1].type = 0;
	nuniforms++;
	// uniforms refering to local variables must be in the 'externals' local list
	assert (g->global == 1 || g->used);
      }
    }
  }
  
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
static Shader * setup_shader (ForeachData * loop, const RegionParameters * region,
			      External * externals,
			      const char * kernel)
{
  for (scalar s in baseblock)
    s.gpu.index = 0;
  int index = 1;
  for (scalar s in baseblock)
    if ((s.input || s.output) && !s.gpu.index) {
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
    s.boundary_left   = s.boundary[left];
    s.boundary_right  = s.boundary[right];
    s.boundary_top    = s.boundary[top];
    s.boundary_bottom = s.boundary[bottom];
  }
  apply_bc_list = loop->dirty;
  
  /**
  ## Allocate reduction fields */

  for (External * g = externals; g && g->name; g++)
    if (g->reduct) {
      scalar s = new scalar;
      s.output = 1;
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
    
  gpu_cpu_sync (baseblock, GL_MAP_WRITE_BIT, loop->fname, loop->line);

  /**
  ## Apply boundary conditions
  
  This can be required if boundary conditions have been modified
  between loops. */

  scalar * listc = NULL;
  for (scalar s in loop->listc)
    if (s.dirty)
      listc = list_prepend (listc, s);
  scalar * listf_x = NULL, * listf_y = NULL;
  foreach_dimension()
    for (scalar s in loop->listf.x)
      if (s.dirty)
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
      s.dirty = false;
      s.gpu.stored = -1;
    }
    free (listc); 
    foreach_dimension() {
      for (scalar s in listf_x) {
	s.dirty = false;
	s.gpu.stored = -1;
      }
      free (listf_x);
    }
    apply_bc_list = loop->dirty;
  }
  
  /**
  For the Intel driver, it looks like the next line is necessary to
  ensure proper synchronisation of the compute shader and fragment
  shader (for example when using output_ppm() for interactive
  display). The nvidia driver somehow does not need this... */

  if (shader->id != GPUContext.current_shader) {
    GL_C (glBindBufferBase (GL_SHADER_STORAGE_BUFFER, 0, 0));
    GL_C (glUseProgram (shader->id));
    for (int i = 0; i < GPUContext.nssbo; i++)
      GL_C (glBindBufferBase (GL_SHADER_STORAGE_BUFFER, i, GPUContext.ssbo[i]));
    GPUContext.current_shader = shader->id;
  }
    
  /**
  ## Set uniforms */

  for (const MyUniform * g = shader->uniforms; g && g->type; g++) {
    void * pointer = g->pointer;
    if (!pointer) {
      assert (g->local >= 0);
      pointer = externals[g->local].pointer;
    }
    switch (g->type) {
    case sym_INT:
      glUniform1iv (g->location, g->nd, pointer); break;
    case sym_FLOAT:
      glUniform1fv (g->location, g->nd, pointer); break;
    case sym_VEC4:
      glUniform4fv (g->location, g->nd, pointer); break;
    case sym_BOOL: {
      int p[g->nd];
      bool * data = pointer;
      for (int i = 0; i < g->nd; i++)
	p[i] = data[i];
      glUniform1iv (g->location, g->nd, p);
      break;
    }
    case sym_LONG: {
      int p[g->nd];
      long * data = pointer;
      for (int i = 0; i < g->nd; i++)
	p[i] = data[i];
      glUniform1iv (g->location, g->nd, p);
      break;
    }
#if SINGLE_PRECISION
    case sym_DOUBLE: {
      float p[g->nd];
      double * data = pointer;
      for (int i = 0; i < g->nd; i++)
	p[i] = data[i];
      glUniform1fv (g->location, g->nd, p);
      break;
    }
    case sym__COORD: {
      float p[2*g->nd];
      double * data = pointer;
      for (int i = 0; i < 2*g->nd; i++)
	p[i] = data[i];
      glUniform2fv (g->location, g->nd, p);
      break;
    }
    case sym_COORD: {
      float p[3*g->nd];
      double * data = pointer;
      for (int i = 0; i < 3*g->nd; i++)
	p[i] = data[i];
      glUniform3fv (g->location, g->nd, p);
      break;
    }
#else // DOUBLE_PRECISION
    case sym_DOUBLE:
      glUniform1dv (g->location, g->nd, pointer); break;
    case sym__COORD:
      glUniform2dv (g->location, g->nd, pointer); break;
    case sym_COORD:
      glUniform3dv (g->location, g->nd, pointer); break;
#endif // DOUBLE_PRECISION
    default:
      assert (false);
    }
  }

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
  ## Render 

  If this is a `foreach_point()` iteration, we draw a single point */

  int Nl = region->level > 0 ? 1 << (region->level - 1) : N/Dimensions.x;
  if (region->n.x == 1 && region->n.y == 1) {
    int csOrigin[] = {
      (region->p.x - X0)/L0*Nl*Dimensions.x,
      (region->p.y - Y0)/L0*Nl*Dimensions.x
    };
    GL_C (glUniform2iv (0, 1, csOrigin));
    assert (!GPUContext.fragment_shader);
    GL_C (glMemoryBarrier (GL_SHADER_STORAGE_BARRIER_BIT));
    GL_C (glDispatchCompute (1, 1, 1));
  }

  /**
  This is a region */
  
  else if (region->n.x || region->n.y) {
    float vsScale[] = {
      (region->box[1].x - region->box[0].x)/L0,
      (region->box[1].y - region->box[0].y)/L0
    };
    float vsOrigin[] = { (region->box[0].x - X0)/L0, (region->box[0].y - Y0)/L0 };
    GL_C (glUniform2fv (1, 1, vsOrigin));
    GL_C (glUniform2fv (2, 1, vsScale));
    assert (GPUContext.fragment_shader);
    GL_C (glMemoryBarrier (GL_SHADER_STORAGE_BARRIER_BIT));
    GL_C (glDrawArrays (GL_TRIANGLES, 0, 6));
  }

  else {
    assert (!GPUContext.fragment_shader);
    GL_C (glMemoryBarrier (GL_SHADER_STORAGE_BARRIER_BIT));
    GL_C (glDispatchCompute (shader->ng[0], shader->ng[1], 1));
  }

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
      s.dirty = false;
    free (loop->dirty), loop->dirty = NULL;
  }
  else {
    gpu_cpu_sync (baseblock, GL_MAP_READ_BIT, loop->fname, loop->line);
    boundary_stencil (loop);
  }
  
  for (scalar s in baseblock)
    if (s.output)
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
