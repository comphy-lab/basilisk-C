/**
# Metal backend for GPUs

This is the [Metal](https://developer.apple.com/metal/) equivalent of
the [CUDA](/src/grid/cuda/cuda.c), [HIP](/src/grid/hip/hip.c) and
[OpenCL](/src/grid/opencl/opencl.c) backends in Basilisk. It is the
only GPU backend which works on recent Apple hardware, since macOS
deprecated OpenGL (version 4.3 and compute shaders were never
supported) and never supported CUDA or a modern OpenCL.

Metal is available on any Mac running macOS >= 10.13 and requires only
the Xcode command-line tools (i.e. no extra library needs to be
installed):

~~~bash
xcode-select --install
cd $BASILISK/grid/metal
make libmetal.a
~~~

The [default Makefile](/src/Makefile.defs) provides recipes to run a
program using the 2D or 3D metal grids:

~~~bash
make myprogram.metal.tst
make myprogram.metal3D.tst
~~~

The test cases are:

~~~bash
cd $BASILISK/test
make metal-tests metal3D-tests
~~~

Note that this file is Objective-C (the Metal API is only available in
Objective-C or C++) and is compiled using `-x objective-c`, see the
[Makefile](Makefile). It is named `metal.c` rather than `metal.m` so
that it fits within the standard Basilisk documentation and
dependency machinery.

## Notes and limitations

* The Metal Shading Language (MSL) does not support double-precision
  floating point numbers. The Metal grids can thus only be used in
  single precision (this is also the default for the other GPU
  backends).
* Fields are stored in a single `MTLBuffer` allocated with
  `MTLResourceStorageModeShared`. On Apple Silicon (where CPU and GPU
  share the same physical memory) this makes the GPU/CPU
  synchronisations in [gpu_cpu_sync_scalar()](#gpu_cpu_sync_scalar) a
  simple `memcpy()` (and no PCIe transfer).
* Unlike the CUDA and OpenCL backends, this backend does not implement
  a compilation cache of its own: `MTLLibrary` objects built from
  source cannot be serialised. None is needed however, since macOS
  caches compiled shaders system-wide (in `MTLCompilerService`), keyed
  on the source. Compiling the ~130 kernels of e.g. [the lid-driven
  cavity](/src/test/lid.c) takes about 3 seconds the first time and
  7 milliseconds on subsequent runs. Kernels are also cached in memory
  for the duration of a run (see
  [load_shader()](/src/grid/gpu/grid.h)).
* `foreach_region()` (i.e. the "fragment shader" mode of the [OpenGL
  backend](/src/grid/gpu/grid.h)) is not implemented, like in the CUDA,
  HIP and OpenCL backends.
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

typedef struct { double x, y, z; } coord;
typedef struct { double x, y, z; } _coord;
typedef struct { float r, g, b, a; } vec4;
typedef float real;
typedef struct { int i; } scalar;
extern int datasize;

typedef struct {
  coord p, * box, n; // region
  int level; // level
} RegionParameters;

extern int N;
extern double X0, Y0, Z0, L0;
extern struct { int x, y, z; } Dimensions;

#include "../../ast/symbols.h"

enum typedef_kind_t {
  sym_SCALAR = sym_root + 1,
  sym_VECTOR,
  sym_TENSOR,
  sym_COORD,
  sym__COORD,
  sym_VEC4,
  sym_IVEC
};

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

#define sysrealloc realloc

#define _GPU 1
#define _CUDA 1
#define _METAL 1
#include "../externals.h"
#include "../gpu/backend.h"

/**
The maximum size of the data which can be passed to a kernel using
`setBytes:length:atIndex:` (larger buffers of uniforms must go through
a `MTLBuffer`). */

#define MAX_SETBYTES 4096

/**
## Metal globals */

static id<MTLDevice> mtl_device = nil;
static id<MTLCommandQueue> mtl_queue = nil;
static id<MTLBuffer> ssbo = nil;

/**
Kernel launches are accumulated within a single command buffer which
is committed (and waited upon) only when the CPU needs to access the
results i.e. in [gpu_synchronize()](#gpu_synchronize). */

static id<MTLCommandBuffer> mtl_cmd = nil;
static id<MTLComputeCommandEncoder> mtl_encoder = nil;

/**
Whether the device supports "non-uniform threadgroups"
i.e. `dispatchThreads:threadsPerThreadgroup:`. This is the case of any
Mac running macOS >= 10.15. */

static bool mtl_nonuniform = false;

#define MTL_CHECK(cond, ...)                                    \
  do {                                                          \
    if (!(cond)) {                                              \
      fprintf (stderr, "%s:%d: Metal error: ", __FILE__, __LINE__); \
      fprintf (stderr, __VA_ARGS__);                            \
      fputc ('\n', stderr);                                     \
      exit (1);                                                 \
    }                                                           \
  } while (0)

typedef struct {
  int type, nd, local;
  size_t size, esize;
  void * pointer;
} MyUniform;

/**
The first field matches the `ulong _data` placeholder of the `_Ctx_`
structure generated for the shader (Metal does not allow pointers in
`constant` buffers). It is only used to detect changes of the field
buffer. */

struct _Ctx_ {
  void * _data;
};

struct _Shader {
  unsigned ng[3], nwg[3];
  int dimension;
  struct _Ctx_ * hctx, * tmp;
  id<MTLBuffer> dctx, dargs;
  MyUniform * uniforms, * locals;
  char * args;
  size_t size, lsize;
  id<MTLLibrary> library;
  id<MTLFunction> function;
  id<MTLComputePipelineState> pipeline;
};

void free_shader (Shader * s)
{
  [s->pipeline release];
  [s->function release];
  [s->library release];
  [s->dctx release];
  [s->dargs release];
  free (s->hctx);
  free (s->tmp);
  free (s->uniforms);
  free (s->locals);
  free (s->args);
  free (s);
}

/**
## Command buffer management */

static id<MTLComputeCommandEncoder> metal_encoder()
{
  if (!mtl_encoder) {
    mtl_cmd = [[mtl_queue commandBuffer] retain];
    MTL_CHECK (mtl_cmd, "could not create command buffer");
    mtl_encoder = [[mtl_cmd computeCommandEncoder] retain];
    MTL_CHECK (mtl_encoder, "could not create compute command encoder");
  }
  return mtl_encoder;
}

static void metal_flush()
{
  if (mtl_encoder) {
    [mtl_encoder endEncoding];
    [mtl_cmd commit];
    [mtl_cmd waitUntilCompleted];
    NSError * err = mtl_cmd.error;
    if (err)
      fprintf (stderr, "%s:%d: Metal warning: %s\n", __FILE__, __LINE__,
               err.localizedDescription.UTF8String);
    [mtl_encoder release], mtl_encoder = nil;
    [mtl_cmd release], mtl_cmd = nil;
  }
}

void gpu_synchronize()
{
  if (mtl_queue)
    @autoreleasepool {
      metal_flush();
    }
}

/**
## Kernel compilation

MSL compilation errors are of the form

~~~
program_source:368:78: error: ...
~~~

which is parsed by [gpu_errors()](/src/grid/gpu/errors.lex) to give
error messages referring to the original Basilisk source code. */

Shader * load_normal_shader (const char * fs,
                             const char * func, const char * file, int line)
{
  if (false) str_append_array (NULL, NULL); // just to avoid a -Wunused-function

  Shader * shader = NULL;
  @autoreleasepool {
    NSString * source = [NSString stringWithUTF8String:fs];
    MTL_CHECK (source, "could not convert the shader source to UTF-8");

    MTLCompileOptions * options = [[MTLCompileOptions alloc] init];
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
    options.fastMathEnabled = YES;
#pragma clang diagnostic pop

    NSError * error = nil;
    id<MTLLibrary> library = [mtl_device newLibraryWithSource:source
                                                     options:options
                                                       error:&error];
    [options release];

    if (!library) {
      if (error) {
        char * merged = gpu_errors (error.localizedDescription.UTF8String, fs,
                                    NULL, "Metal");
        fputs (merged, stderr);
        free (merged);
      }
      return NULL;
    }

    id<MTLFunction> function = [library newFunctionWithName:
                                          [NSString stringWithUTF8String:func]];
    if (!function) {
      fprintf (stderr, "%s:%d: Metal: error: could not find kernel '%s'\n",
               file, line, func);
      [library release];
      return NULL;
    }

    error = nil;
    id<MTLComputePipelineState> pipeline =
      [mtl_device newComputePipelineStateWithFunction:function error:&error];
    if (!pipeline) {
      fprintf (stderr, "%s:%d: Metal: error: %s\n", file, line,
               error ? error.localizedDescription.UTF8String : "unknown error");
      [function release];
      [library release];
      return NULL;
    }

    shader = calloc (1, sizeof (Shader));
    shader->library = library;
    shader->function = function;
    shader->pipeline = pipeline;
  }
  return shader;
}

/**
## Context and memory management */

bool gpu_init_context (GPUData ** data)
{
  bool initialized = (mtl_device != nil);
  if (!initialized) {
    mtl_device = MTLCreateSystemDefaultDevice();
    MTL_CHECK (mtl_device, "no Metal device found");
    mtl_queue = [mtl_device newCommandQueue];
    MTL_CHECK (mtl_queue, "could not create a Metal command queue");
    mtl_nonuniform = ([mtl_device supportsFamily: MTLGPUFamilyApple4] ||
                      [mtl_device supportsFamily: MTLGPUFamilyMac2]);
  }
  *data = NULL;
  return !initialized;
}

void gpu_free_context (GPUData * data)
{
  gpu_synchronize();
  [ssbo release], ssbo = nil;
  GPUContext.current_size = 0;
}

void realloc_ssbo (size_t field_size)
{
  if (!datasize)
    return;
  size_t totalsize = field_size*datasize;
  assert (totalsize > GPUContext.current_size);
  gpu_synchronize();
  if (ssbo)
    [ssbo release];
  ssbo = [mtl_device newBufferWithLength:totalsize
                     options:MTLResourceStorageModeShared];
  MTL_CHECK (ssbo, "could not allocate %ld bytes on the GPU", totalsize);
  GPUContext.current_size = totalsize;
}

/**
Since the field buffer uses `MTLResourceStorageModeShared`, GPU/CPU
synchronisations are just memory copies (and are free on systems with
a unified memory architecture). */

void gpu_cpu_sync_scalar (int i, int block, char * data, size_t field_size, SyncMode mode)
{
  size_t size = field_size*sizeof(real), offset = i*size, totalsize = block*size;
  gpu_synchronize();
  char * cd = data + offset, * gd = (char *) ssbo.contents + offset;
  if (mode == GPU_READ)
    memcpy (cd, gd, totalsize);
  else if (mode == GPU_WRITE)
    memcpy (gd, cd, totalsize);
  else
    assert (false);
}

void reset_scalar (int i, int block, size_t field_size, double val)
{
  size_t offset = i*field_size, n = max(block, 1)*field_size;
  gpu_synchronize();
  real * p = (real *) ssbo.contents + offset, v = val;
  if (!v)
    memset (p, 0, n*sizeof (real));
  else
    for (size_t j = 0; j < n; j++)
      p[j] = v;
}

/**
## Uniforms

This is identical to the [CUDA backend](/src/grid/cuda/cuda.c): global
variables are gathered within a single `_Ctx_` structure and local
variables within a `_Local` structure, both of which are passed as
kernel arguments. */

// Calculate padding needed to align 'current_offset' to 'alignment'
static size_t pad_to_align (size_t current_offset, size_t alignment) {
  return (alignment - (current_offset % alignment)) % alignment;
}

void finalize_shader (Shader * shader, External * externals, External * merged,
                      unsigned ng[3], unsigned nwg[3], int dim)
{
  for (int i = 0; i < 3; i++)
    shader->ng[i] = ng[i], shader->nwg[i] = nwg[i];
  shader->dimension = dim;

  /**
  ## Make list of local and global uniforms */

  for (External * g = merged; g; g = g->next)
    g->used = 0;
  int index = 1;
  for (External * g = externals; g && g->name; g++)
    g->used = index++;
  int nglobals = 0, nlocals = 0;
  for (const External * g = merged; g; g = g->next)
    if (is_external_variable (g)) {
      int nd = g->data ? ((int *)g->data)[0] : 1;
      size_t esize;
      switch (g->type) {
      case sym_INT: case sym_LONG:
        esize = sizeof(int); break;
      case sym_FLOAT:
        esize = sizeof(float); break;
      case sym_VEC4:
        nd *= 4;
        esize = sizeof(float); break;
      case sym_BOOL:
        esize = sizeof(bool); break;

      /**
      Metal does not support double-precision floating point numbers:
      `real` is always `float`. */

      case sym_DOUBLE:
        esize = sizeof(float); break;
      case sym__COORD:
        nd *= dim;
        esize = sizeof(float); break;
      case sym_COORD:
        nd *= 3;
        esize = sizeof(float); break;
      default:
        assert (false);
      }
      MyUniform uniform = {
        .type = g->type, .nd = nd, .size = nd*esize, .esize = esize,
        .local = g->global == 1 ? -1 : g->used - 1,
        .pointer = g->global == 1 ? g->pointer : NULL
      };
      if (g->global) {
        shader->uniforms = realloc (shader->uniforms, (nglobals + 2)*sizeof(MyUniform));
        shader->uniforms[nglobals] = uniform;
        shader->uniforms[nglobals + 1].type = 0;
        nglobals++;
      }
      else {
        shader->locals = realloc (shader->locals, (nlocals + 2)*sizeof(MyUniform));
        shader->locals[nlocals] = uniform;
        shader->locals[nlocals + 1].type = 0;
        nlocals++;
      }
    }

  /**
  Allocate the host buffers holding the uniforms. A device buffer is
  necessary only for the (unlikely) case where they do not fit within
  the `setBytes:length:atIndex:` limit. */

  shader->size = sizeof (struct _Ctx_);
  for (const MyUniform * g = shader->uniforms; g && g->type; g++)
    shader->size += pad_to_align (shader->size, g->esize) + g->size;
  shader->hctx = calloc (1, shader->size);
  shader->tmp = calloc (1, shader->size);
  if (shader->size > MAX_SETBYTES) {
    shader->dctx = [mtl_device newBufferWithLength:shader->size
                               options:MTLResourceStorageModeShared];
    MTL_CHECK (shader->dctx, "could not allocate %ld bytes on the GPU", shader->size);
  }

  /**
  Allocate host buffer to hold locals. */

  if (shader->locals) {
    for (const MyUniform * g = shader->locals; g && g->type; g++)
      shader->lsize += pad_to_align (shader->lsize, g->esize) + g->size;
    shader->args = calloc (1, shader->lsize);
    if (shader->lsize > MAX_SETBYTES) {
      shader->dargs = [mtl_device newBufferWithLength:shader->lsize
                                  options:MTLResourceStorageModeShared];
      MTL_CHECK (shader->dargs, "could not allocate %ld bytes on the GPU", shader->lsize);
    }
  }
}

static char * set_uniforms (const MyUniform * uniforms,
                            const External * externals,
                            char * buffer, char * start)
{
  for (const MyUniform * g = uniforms; g && g->type; g++) {
    void * pointer = g->pointer;
    if (!pointer) {
      assert (g->local >= 0);
      pointer = externals[g->local].pointer;
    }
    buffer += pad_to_align (buffer - start, g->esize);
    switch (g->type) {
    case sym_INT: case sym_FLOAT: case sym_VEC4: case sym_BOOL:
      memcpy (buffer, pointer, g->size);
      break;
    case sym_LONG: {
      int p[g->nd];
      long * data = pointer;
      for (int i = 0; i < g->nd; i++)
	p[i] = data[i];
      memcpy (buffer, p, g->size);
      break;
    }
    case sym_DOUBLE: case sym__COORD: case sym_COORD: {
      float p[g->nd];
      double * data = pointer;
      for (int i = 0; i < g->nd; i++)
	p[i] = data[i];
      memcpy (buffer, p, g->size);
      break;
    }
    default:
      assert (false);
    }
    buffer += g->size;
  }
  return buffer;
}

void post_setup_shader (Shader * shader, External * externals)
{
  /**
  Set the (placeholder for the) field buffer. */

  assert (ssbo);
  shader->tmp->_data = (__bridge void *) ssbo;
  char * buffer = ((char *)shader->tmp) + sizeof(struct _Ctx_);

  /**
  Set globals */

  buffer = set_uniforms (shader->uniforms, externals, buffer, (char *)shader->tmp);
  assert (shader->size == buffer - (char *)shader->tmp);

  if (memcmp (shader->hctx, shader->tmp, shader->size)) {
    struct _Ctx_ * tmp = shader->tmp; shader->tmp = shader->hctx; shader->hctx = tmp;
    if (shader->dctx) {
      gpu_synchronize(); // the buffer may still be in use by a pending kernel
      memcpy (shader->dctx.contents, shader->hctx, shader->size);
    }
  }

  /**
  Set locals */

  buffer = set_uniforms (shader->locals, externals, shader->args, shader->args);
  assert (shader->lsize == buffer - shader->args);
  if (shader->dargs) {
    gpu_synchronize();
    memcpy (shader->dargs.contents, shader->args, shader->lsize);
  }
}

/**
## Kernel launch */

static void metal_dispatch (id<MTLComputeCommandEncoder> encoder,
                            id<MTLComputePipelineState> pipeline,
                            NSUInteger gx, NSUInteger gy, NSUInteger gz,
                            NSUInteger tx, NSUInteger ty, NSUInteger tz)
{
  /**
  Register-hungry kernels may not be able to use the default
  threadgroup size. */

  NSUInteger maxt = pipeline.maxTotalThreadsPerThreadgroup;
  while (tx*ty*tz > maxt && (tx > 1 || ty > 1 || tz > 1)) {
    if (tz >= ty && tz >= tx) tz = (tz + 1)/2;
    else if (ty >= tx) ty = (ty + 1)/2;
    else tx = (tx + 1)/2;
  }
  if (mtl_nonuniform)
    [encoder dispatchThreads: MTLSizeMake (gx, gy, gz)
             threadsPerThreadgroup: MTLSizeMake (tx, ty, tz)];
  else
    [encoder dispatchThreadgroups: MTLSizeMake ((gx + tx - 1)/tx, (gy + ty - 1)/ty,
                                                (gz + tz - 1)/tz)
             threadsPerThreadgroup: MTLSizeMake (tx, ty, tz)];
}

int run_shader (const Shader * shader, const RegionParameters * region)
{
  int Nl = region->level > 0 ? 1 << (region->level - 1) : N/Dimensions.x;

  @autoreleasepool {
    id<MTLComputeCommandEncoder> encoder = metal_encoder();
    [encoder setComputePipelineState: shader->pipeline];
    [encoder setBuffer: ssbo offset: 0 atIndex: 0];
    if (shader->dctx)
      [encoder setBuffer: shader->dctx offset: 0 atIndex: 1];
    else
      [encoder setBytes: shader->hctx length: shader->size atIndex: 1];
    if (shader->args) {
      if (shader->dargs)
        [encoder setBuffer: shader->dargs offset: 0 atIndex: 3];
      else
        [encoder setBytes: shader->args length: shader->lsize atIndex: 3];
    }

    struct { int x, y, z; } csOrigin = {0,0,0};
    size_t csOriginSize = shader->dimension*sizeof(int);

    /**
    ## Render

    If this is a `foreach_point()` iteration, we access a single point */

    if (region->n.x == 1 && region->n.y == 1 &&
        (shader->dimension == 2 || region->n.z == 1)) {
      assert (!GPUContext.fragment_shader);
      csOrigin.x = (region->p.x - X0)/L0*Nl*Dimensions.x;
      csOrigin.y = (region->p.y - Y0)/L0*Nl*Dimensions.x;
      if (shader->dimension == 3)
        csOrigin.z = (region->p.z - Z0)/L0*Nl*Dimensions.x;
      [encoder setBytes: &csOrigin length: csOriginSize atIndex: 2];
      metal_dispatch (encoder, shader->pipeline, 1, 1, 1, 1, 1, 1);
    }

    /**
    This is a region */

    else if (region->n.x || region->n.y ||
             (shader->dimension == 3 && region->n.z)) {

      /**
      A region requires the "fragment shader" mode of the [OpenGL
      backend](/src/grid/gpu/grid.h), which has no Metal equivalent
      (this is also the case of the CUDA, HIP and OpenCL backends). */

      fprintf (stderr, "%s:%d: Metal: error: foreach_region() is not implemented\n",
               __FILE__, __LINE__);
      exit (1);
    }

    else {
      assert (!GPUContext.fragment_shader);
      [encoder setBytes: &csOrigin length: csOriginSize atIndex: 2];
      metal_dispatch (encoder, shader->pipeline,
                      shader->ng[0]*shader->nwg[0], shader->ng[1]*shader->nwg[1],
                      shader->ng[2]*shader->nwg[2],
                      shader->nwg[0], shader->nwg[1], shader->nwg[2]);
    }
  }
  return Nl;
}

void gpu_free_solver (void)
{
  gpu_synchronize();
  [mtl_queue release], mtl_queue = nil;
  mtl_device = nil; // owned by the system
}

/**
## Reductions

A standard two-elements-per-thread, multi-pass tree reduction, as in
the [CUDA](/src/grid/cuda/cuda.c) and
[OpenCL](/src/grid/opencl/opencl.c) backends. */

#define REDUCE_THREADS 256

static const char reduce_source[] =
  "#include <metal_stdlib>\n"
  "using namespace metal;\n"
  "#define REDUCE(reduced,rhs) %s\n"
  "kernel void reduce (device const float * input [[buffer(0)]],\n"
  "                    device float * output [[buffer(1)]],\n"
  "                    constant uint & n [[buffer(2)]],\n"
  "                    constant uint & offset [[buffer(3)]],\n"
  "                    uint tid [[thread_position_in_threadgroup]],\n"
  "                    uint gid [[threadgroup_position_in_grid]],\n"
  "                    uint nt [[threads_per_threadgroup]]) {\n"
  "  threadgroup float sdata[%d];\n"
  "  uint i = gid*nt*2 + tid;\n"
  "  float reduced = %s;\n"
  "  if (i < n)\n"
  "    REDUCE (reduced, input[i + offset]);\n"
  "  if (i + nt < n)\n"
  "    REDUCE (reduced, input[i + nt + offset]);\n"
  "  sdata[tid] = reduced;\n"
  "  threadgroup_barrier (mem_flags::mem_threadgroup);\n"
  "  for (uint s = nt/2; s > 0; s >>= 1) {\n"
  "    if (tid < s)\n"
  "      REDUCE (sdata[tid], sdata[tid + s]);\n"
  "    threadgroup_barrier (mem_flags::mem_threadgroup);\n"
  "  }\n"
  "  if (tid == 0)\n"
  "    output[gid] = sdata[0];\n"
  "}\n";

static id<MTLComputePipelineState> compile_kernel (const char * start, const char * op)
{
  id<MTLComputePipelineState> pipeline = nil;
  @autoreleasepool {
    char source[sizeof (reduce_source) + 200];
    snprintf (source, sizeof (source), reduce_source, op, REDUCE_THREADS, start);
    NSError * error = nil;
    id<MTLLibrary> library =
      [mtl_device newLibraryWithSource: [NSString stringWithUTF8String: source]
                  options: nil error: &error];
    MTL_CHECK (library, "could not compile the reduction kernel: %s",
               error ? error.localizedDescription.UTF8String : "unknown error");
    id<MTLFunction> function = [library newFunctionWithName: @"reduce"];
    MTL_CHECK (function, "could not find the reduction kernel");
    error = nil;
    pipeline = [mtl_device newComputePipelineStateWithFunction: function error: &error];
    MTL_CHECK (pipeline, "could not create the reduction pipeline: %s",
               error ? error.localizedDescription.UTF8String : "unknown error");
    [function release];
    [library release];
  }
  return pipeline;
}

static float metal_reduce (id<MTLBuffer> d_input, const size_t nb, const char op,
                           const size_t initial_offset)
{
  id<MTLComputePipelineState> pipeline;
  switch (op) {
  case '+': {
    static id<MTLComputePipelineState> k = nil;
    if (!k) k = compile_kernel ("0.0f", "reduced += rhs");
    pipeline = k;
    break;
  }
  case 'm': {
    static id<MTLComputePipelineState> k = nil;
    if (!k) k = compile_kernel ("3e38f", "reduced = min(reduced,rhs)");
    pipeline = k;
    break;
  }
  case 'M': {
    static id<MTLComputePipelineState> k = nil;
    if (!k) k = compile_kernel ("-3e38f", "reduced = max(reduced,rhs)");
    pipeline = k;
    break;
  }
  default:
    assert (0);
    return 0.;
  }

  /**
  The scratch buffers are reallocated only when necessary. */

  static size_t Np = 0;
  static id<MTLBuffer> d_output_a = nil, d_output_b = nil;
  if (nb > Np) {
    const size_t max_blocks = (nb + 2*REDUCE_THREADS - 1)/(2*REDUCE_THREADS) + 1;
    [d_output_a release];
    [d_output_b release];
    d_output_a = [mtl_device newBufferWithLength: max_blocks*sizeof(float)
                             options: MTLResourceStorageModeShared];
    d_output_b = [mtl_device newBufferWithLength: max_blocks*sizeof(float)
                             options: MTLResourceStorageModeShared];
    MTL_CHECK (d_output_a && d_output_b, "could not allocate the reduction buffers");
    Np = nb;
  }

  /**
  Multi-pass reduction. Only the first pass applies the offset within
  the field buffer, subsequent passes read tightly-packed
  intermediate results. */

  @autoreleasepool {
    id<MTLBuffer> input = d_input, output = d_output_a;
    size_t current_n = nb, current_offset = initial_offset;
    id<MTLComputeCommandEncoder> encoder = metal_encoder();
    [encoder setComputePipelineState: pipeline];
    while (current_n > 1) {
      uint32_t blocks = (current_n + 2*REDUCE_THREADS - 1)/(2*REDUCE_THREADS);
      uint32_t n_param = current_n, offset_param = current_offset;
      [encoder setBuffer: input offset: 0 atIndex: 0];
      [encoder setBuffer: output offset: 0 atIndex: 1];
      [encoder setBytes: &n_param length: sizeof (n_param) atIndex: 2];
      [encoder setBytes: &offset_param length: sizeof (offset_param) atIndex: 3];
      [encoder dispatchThreadgroups: MTLSizeMake (blocks, 1, 1)
               threadsPerThreadgroup: MTLSizeMake (REDUCE_THREADS, 1, 1)];
      current_n = blocks;
      input = output;
      output = (output == d_output_a) ? d_output_b : d_output_a;
      current_offset = 0;
    }
    metal_flush();
    return *((float *) input.contents);
  }
}

double gpu_reduction (size_t offset,
                      const char op,
                      const RegionParameters * region,
                      GPUData * data,
                      size_t nb,
                      int dim)
{
  if (region->n.x == 1 && region->n.y == 1 && (dim == 2 || region->n.z == 1)) {
    int i = (region->p.x - X0)/L0*N;
    int j = (region->p.y - Y0)/L0*N;
    int k = dim == 3 ? (region->p.z - Z0)/L0*N : 0;
    if (i < 0 || i >= N || j < 0 || j >= N || (dim == 3 && (k < 0 || k >= N)))
      return 0.;
    offset += dim == 2 ? i*N + j : (i*N + j)*N + k;
    nb = 1;
  }
  if (nb == 1) {
    gpu_synchronize();
    return ((real *) ssbo.contents)[offset];
  }
  return metal_reduce (ssbo, nb, op, offset);
}

/**
## Device characteristics */

void printWorkGroupsCapabilities()
{
  gpu_limits (stderr);
}

void gpu_limits (FILE * fp)
{
  if (!mtl_device) {
    fprintf (fp, "Metal: no device\n");
    return;
  }
  @autoreleasepool {
    fprintf (fp, "Metal device:                    %s\n",
             mtl_device.name.UTF8String);
    fprintf (fp, "Unified memory:                  %s\n",
             mtl_device.hasUnifiedMemory ? "yes" : "no");
    fprintf (fp, "Non-uniform threadgroups:        %s\n",
             mtl_nonuniform ? "yes" : "no");
    fprintf (fp, "Recommended working set size:    %llu\n",
             (unsigned long long) mtl_device.recommendedMaxWorkingSetSize);
    fprintf (fp, "Maximum buffer length:           %lu\n",
             (unsigned long) mtl_device.maxBufferLength);
    fprintf (fp, "Maximum threads per threadgroup: %lu x %lu x %lu\n",
             (unsigned long) mtl_device.maxThreadsPerThreadgroup.width,
             (unsigned long) mtl_device.maxThreadsPerThreadgroup.height,
             (unsigned long) mtl_device.maxThreadsPerThreadgroup.depth);
    fprintf (fp, "Maximum threadgroup memory:      %lu\n",
             (unsigned long) mtl_device.maxThreadgroupMemoryLength);
  }
}

/**
## See also

* [Metal Shading Language Specification](https://developer.apple.com/metal/Metal-Shading-Language-Specification.pdf)
* [Grids on GPUs](/src/grid/gpu/grid.h)
* [The CUDA backend](/src/grid/cuda/cuda.c)
* [The OpenCL backend](/src/grid/opencl/opencl.c)
*/
