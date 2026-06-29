/**
# OpenCL backend for GPUs
 
This is the [OpenCL](https://www.khronos.org/opencl/) equivalent of
the [CUDA backend](/src/grid/cuda/cuda.c) in Basilisk.  It uses the
OpenCL runtime API, which can be installed on Debian systems using:

~~~bashrc
apt install ocl-icd-opencl-dev clinfo
cd $BASILISK/grid/opencl
make libocl.a
~~~

The [default Makefile](/src/Makefile.defs) provides a recipe to run a
program using the opencl grid:

~~~bashrc
make myprogram.ocl.tst
~~~

The tests cases are:

~~~bashrc
cd $BASILISK/test
make ocl-tests
~~~
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <sys/stat.h>
#include <errno.h>
#include "a32.h"

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
extern struct { int x, y; } Dimensions;

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
#define _OPENCL 1
#include "../externals.h"
#include "../gpu/backend.h"

#define CL_TARGET_OPENCL_VERSION 200
#include <CL/cl.h>

// OpenCL globals
static cl_mem ssbo = NULL;
static cl_device_id dev = NULL;
static cl_context ctx = NULL;
static cl_command_queue queue = NULL;
static cl_platform_id platform = NULL;

#define CL_CHECK(x) \
  do { \
    cl_int _err = (x); \
    if (_err != CL_SUCCESS) { \
      fprintf(stderr, "%s:%d: OpenCL error: %d\n", __FILE__, __LINE__, _err); \
      exit(1); \
    } \
  } while(0)

typedef struct {
  int type, nd, local;
  size_t size, esize;
  void * pointer;
} MyUniform;

struct _Ctx_ {
  cl_mem _data;
};

struct _Shader {
  unsigned ng[2], nwg[2];
  struct _Ctx_ * hctx, * tmp;
  cl_mem dctx;
  MyUniform * uniforms, * locals;
  char * args;
  size_t size, lsize;
  cl_kernel kernel;
  cl_program program;
};

void free_shader (Shader * s)
{
  if (s->kernel) clReleaseKernel(s->kernel);
  if (s->program) clReleaseProgram(s->program);
  if (s->dctx) clReleaseMemObject(s->dctx);
  free (s->hctx);
  free (s->tmp);
  free (s->uniforms);
  free (s->locals);
  free (s->args);
  free (s);
}

static
void architecture (char * arch)
{
  // OpenCL doesn't use sm_XX like CUDA, but we can use generic flags
  sprintf(arch, "-cl-fast-relaxed-math");
}

static
int create_tmpdir (const char * path)
{
  struct stat st;
  if (stat (path, &st) == 0) {
    if (S_ISDIR (st.st_mode))
      return 0; // Directory exists
    else {
      errno = ENOTDIR;
      return -1; // Path exists but is not a directory
    }
  }
  // Directory does not exist, try to create it
  if (mkdir (path, 0755) == 0)
    return 0; // Successfully created
  return -1; // Failed to create
}

Shader * load_normal_shader (const char * fs,
                             const char * func, const char * file, int line)
{
  // ---------------------------------------------------------------
  // Try to read from a compilation cache (by default in /tmp/buda/)
  // ---------------------------------------------------------------
  
  Adler32Hash hasha;
  a32_hash_init (&hasha);
  a32_hash_add (&hasha, fs, strlen (fs));
  uint32_t hash = a32_hash (&hasha);

  const char * tmpdir = getenv ("TMPDIR"), * tmp = tmpdir ? tmpdir : "/tmp";
  char cache[strlen(tmp) + strlen("/buda/ffffffff") + 1];
  sprintf (cache, "%s/buda", tmp);
  if (create_tmpdir (cache)) {
    fprintf (stderr, "%s:%d: cannot create temporary directory '%s'\n", \
             __FILE__, __LINE__, cache);
    perror ("");
    exit (1);
  }
  sprintf (cache, "%s/buda/%x", tmp, hash);
  
  // Check if cached binary exists
  struct stat st;
  cl_program program = NULL;
  
  if (stat (cache, &st) == 0) { // found in cache
    FILE * fp = fopen (cache, "rb");
    if (fp) {
      size_t binary_size = st.st_size;
      unsigned char * binary = malloc (binary_size);
      if (binary && fread (binary, 1, binary_size, fp) == binary_size) {
        cl_int err;
        program = clCreateProgramWithBinary (ctx, 1, &dev, &binary_size, 
                                              (const unsigned char **)&binary, NULL, &err);
        if (err == CL_SUCCESS)
          err = clBuildProgram (program, 1, &dev, "", NULL, NULL);
        if (err != CL_SUCCESS) {
          clReleaseProgram (program);
          program = NULL;
        }
      }
      if (binary) free (binary);
      fclose (fp);
    }
  }
  
  if (!program) { // not found in cache or failed to load
    char arch[256];
    architecture (arch);
    const char * sources[] = { fs };
    cl_int err;
    program = clCreateProgramWithSource (ctx, 1, sources, NULL, &err);
    CL_CHECK (err);
    err = clBuildProgram (program, 1, &dev, arch, NULL, NULL);
    if (err != CL_SUCCESS) {
      size_t log_size;
      clGetProgramBuildInfo (program, dev, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
      if (log_size > 1) {
        char * log = (char *) malloc (log_size);
        clGetProgramBuildInfo (program, dev, CL_PROGRAM_BUILD_LOG, log_size, log, NULL);
        char * error = gpu_errors (log, fs, NULL, "OpenCL");
        fputs (error, stderr);
        free (error);
        free (log);
      }
      clReleaseProgram (program);
      return NULL;
    }

#if 0
    {
      size_t log_size;
      clGetProgramBuildInfo(program, dev,
                            CL_PROGRAM_BUILD_LOG,
                            0, NULL, &log_size);
      
      char *log = malloc(log_size);
      clGetProgramBuildInfo(program, dev,
                            CL_PROGRAM_BUILD_LOG,
                            log_size, log, NULL);
      
      fprintf (stderr, "%s:%d: %s\n", __FILE__, __LINE__, log);
      free (log);
    }
#endif
    
    // Cache the program binary
    FILE * fp = fopen (cache, "wb");
    if (fp) {
      size_t * binary_sizes = malloc (sizeof (size_t));
      unsigned char ** binaries = malloc (sizeof (unsigned char *));
      if (binary_sizes && binaries) {
        CL_CHECK (clGetProgramInfo (program, CL_PROGRAM_BINARY_SIZES, sizeof (size_t), binary_sizes, NULL));
        binaries[0] = malloc (binary_sizes[0]);
        if (binaries[0])
          CL_CHECK (clGetProgramInfo (program, CL_PROGRAM_BINARIES, sizeof (unsigned char *), binaries, NULL));
        if (binaries[0] && binary_sizes[0] > 0)
          fwrite (binaries[0], 1, binary_sizes[0], fp);
        free (binaries[0]);
      }
      free (binaries);
      free (binary_sizes);
      fclose (fp);
    }
  }
  
  // Create kernel
  cl_int err;
  cl_kernel kernel = clCreateKernel (program, func, &err);
  CL_CHECK (err);

  Shader * shader = calloc (1, sizeof (Shader));
  shader->program = program;
  shader->kernel = kernel;
  return shader;
}

bool gpu_init_context (GPUData ** data)
{
  bool initialized = (ctx != NULL);
  if (!initialized) {
    cl_uint num_platforms, num_devices;
    cl_int err;
    
    // Get platform
    err = clGetPlatformIDs (1, &platform, &num_platforms);
    CL_CHECK (err);
    
    // Get GPU device
    err = clGetDeviceIDs (platform, CL_DEVICE_TYPE_GPU, 1, &dev, &num_devices);
    if (err != CL_SUCCESS) {
      // Try default device
      err = clGetDeviceIDs (platform, CL_DEVICE_TYPE_DEFAULT, 1, &dev, &num_devices);
      CL_CHECK (err);
    }
    
    // Create context
    ctx = clCreateContext (NULL, 1, &dev, NULL, NULL, &err);
    CL_CHECK (err);
    
    // Create command queue with properties for OpenCL 2.0+
    cl_queue_properties properties[] = { CL_QUEUE_PROPERTIES, CL_QUEUE_PROFILING_ENABLE, 0 };
    queue = clCreateCommandQueueWithProperties (ctx, dev, properties, &err);
    CL_CHECK (err);
  }
  *data = NULL;
  return !initialized;
  if (false) str_append_array (NULL, NULL); // just to avoid a -Wunused-function
}

void gpu_free_context (GPUData * data)
{
  if (ssbo) {
    CL_CHECK (clReleaseMemObject (ssbo));
    ssbo = NULL;
  }
  GPUContext.current_size = 0;
}

void realloc_ssbo (size_t field_size)
{
  if (!datasize)
    return;
  size_t totalsize = field_size * datasize;
  assert (totalsize > GPUContext.current_size);
  cl_mem new_buffer = clCreateBuffer (ctx, CL_MEM_READ_WRITE, totalsize, NULL, NULL);
  assert (new_buffer);
  if (GPUContext.current_size > 0) {
    // Copy old data to new buffer
    CL_CHECK (clEnqueueCopyBuffer (queue, ssbo, new_buffer, 0, 0, GPUContext.current_size, 0, NULL, NULL));
    CL_CHECK (clFinish (queue));
    CL_CHECK (clReleaseMemObject (ssbo));
  }
  ssbo = new_buffer;
  GPUContext.current_size = totalsize;
}

void gpu_cpu_sync_scalar (int i, int block, char * data, size_t field_size, SyncMode mode)
{
  size_t size = field_size * sizeof (real);
  size_t offset = i * size;
  size_t totalsize = block * size;
  char * cd = data + offset;
  
  if (mode == GPU_READ) {
    CL_CHECK (clEnqueueReadBuffer (queue, ssbo, CL_TRUE, offset, totalsize, cd, 0, NULL, NULL));
    CL_CHECK (clFinish(queue));    
  }
  else if (mode == GPU_WRITE)
    CL_CHECK (clEnqueueWriteBuffer (queue, ssbo, CL_TRUE, offset, totalsize, cd, 0, NULL, NULL));
  else
    assert (false);
}

void reset_scalar (int i, int block, size_t field_size, double val)
{
  size_t size = field_size * sizeof (real);
  size_t offset = i * size;
  size_t totalsize = max (block, 1) * size;

#if SINGLE_PRECISION
  float fval = val;
  CL_CHECK (clEnqueueFillBuffer (queue, ssbo, &fval, sizeof (float), offset, totalsize, 0, NULL, NULL));
  CL_CHECK (clFinish (queue));
#else
  fprintf (stderr, "%s:%d: error: not implemented yet\n", __FILE__, __LINE__);
#endif
}

// Calculate padding needed to align 'current_offset' to 'alignment'
static size_t pad_to_align (size_t current_offset, size_t alignment) {
  return (alignment - (current_offset % alignment)) % alignment;
}

void finalize_shader (Shader * shader, External * externals, External * merged,
                      unsigned ng[2], unsigned nwg[2])
{
  for (int i = 0; i < 2; i++)
    shader->ng[i] = ng[i], shader->nwg[i] = nwg[i];

  /**
   * Make list of local and global uniforms
   */

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
        esize = sizeof (int); break;
      case sym_FLOAT:
        esize = sizeof (float); break;
      case sym_VEC4:
        nd *= 4;
        esize = sizeof (float); break;
      case sym_BOOL:
        esize = sizeof (int); break;
#if SINGLE_PRECISION
      case sym_DOUBLE:
        esize = sizeof (float); break;
      case sym__COORD:
        nd *= 2;
        esize = sizeof (float); break;
      case sym_COORD:
        nd *= 3;
        esize = sizeof (float); break;
#else // DOUBLE_PRECISION
      case sym_DOUBLE:
        esize = sizeof (double); break;
      case sym__COORD:
        nd *= 2;
        esize = sizeof (double); break;
      case sym_COORD:
        nd *= 3;
        esize = sizeof (double); break;
#endif // DOUBLE_PRECISION
      default:
        assert (false);
      }
      MyUniform uniform = {
        .type = g->type, .nd = nd, .size = nd*esize, .esize = esize,
        .local = g->global == 1 ? -1 : g->used - 1,
        .pointer = g->global == 1 ? g->pointer : NULL
      };
      if (g->global) {
        shader->uniforms = realloc (shader->uniforms, (nglobals + 2)*sizeof (MyUniform));
        shader->uniforms[nglobals] = uniform;
        shader->uniforms[nglobals + 1].type = 0;
        nglobals++;
      }
      else {
        shader->locals = realloc (shader->locals, (nlocals + 2)*sizeof (MyUniform));
        shader->locals[nlocals] = uniform;
        shader->locals[nlocals + 1].type = 0;
        nlocals++;
      }
    }

  /**
   * Allocate host and device buffers to hold uniforms.
   */

  shader->size = sizeof (struct _Ctx_);
  for (const MyUniform * g = shader->uniforms; g && g->type; g++)
    shader->size += pad_to_align (shader->size, g->esize) + g->size;
  shader->hctx = calloc (1, shader->size);
  shader->tmp = calloc (1, shader->size);
  shader->dctx = clCreateBuffer (ctx, CL_MEM_READ_ONLY, shader->size, NULL, NULL);

  /**
   * Allocate host buffer to hold locals.
   */

  if (shader->locals) {
    for (const MyUniform * g = shader->locals; g && g->type; g++)
      shader->lsize += pad_to_align (shader->lsize, g->esize) + g->size;
    shader->args = calloc (1, shader->lsize);
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
    case sym_INT: case sym_FLOAT: case sym_VEC4:
      memcpy (buffer, pointer, g->size);
      break;
    case sym_BOOL: {
      int b = *((bool *)pointer);
      memcpy (buffer, &b, g->size);
      break;
    }
    case sym_LONG: {
      int p[g->nd];
      long * data = pointer;
      for (int i = 0; i < g->nd; i++)
        p[i] = data[i];
      memcpy (buffer, p, g->size);
      break;
    }
#if SINGLE_PRECISION
    case sym_DOUBLE: case sym__COORD: case sym_COORD: {
      float p[g->nd];
      double * data = pointer;
      for (int i = 0; i < g->nd; i++)
        p[i] = data[i];
      memcpy (buffer, p, g->size);
      break;
    }
#else // DOUBLE_PRECISION
    case sym_DOUBLE: case sym__COORD: case sym_COORD:
      memcpy (buffer, pointer, g->size);
      break;
#endif // DOUBLE_PRECISION
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
   * Set SSBO pointer.
   */

  assert (ssbo);
  shader->tmp->_data = ssbo;
  char * buffer = ((char *)shader->tmp) + sizeof (struct _Ctx_);

  /**
   * Set globals
   */

  buffer = set_uniforms (shader->uniforms, externals, buffer, (char *)shader->tmp);
  assert (shader->size == buffer - (char *)shader->tmp);

  if (memcmp (shader->hctx, shader->tmp, shader->size)) {
    struct _Ctx_ * tmp = shader->tmp; shader->tmp = shader->hctx; shader->hctx = tmp;
    CL_CHECK (clEnqueueWriteBuffer (queue, shader->dctx, CL_TRUE, 0, shader->size, shader->hctx, 0, NULL, NULL));
  }

  /**
   * Set locals
   */

  buffer = set_uniforms (shader->locals, externals, shader->args, shader->args);
  assert (shader->lsize == buffer - shader->args);
}

int run_shader (const Shader * shader, const RegionParameters * region)
{
  struct { int x, y; } csOrigin = {0,0};

  /**
   * Set kernel arguments
   */

  CL_CHECK (clFinish (queue)); // fixme: which synchro are necessary?

  // fixme: possible optimisation, only set args which have changed
  CL_CHECK (clSetKernelArg (shader->kernel, 0, sizeof (cl_mem), &ssbo));
  CL_CHECK (clSetKernelArg (shader->kernel, 1, sizeof (cl_mem), &shader->dctx));
  CL_CHECK (clSetKernelArg (shader->kernel, 2, sizeof (struct { int x, y; }), &csOrigin));
  if (shader->args)
    CL_CHECK (clSetKernelArg (shader->kernel, 3, shader->lsize, shader->args));

  /**
   * Render
   * If this is a foreach_point() iteration, we access a single point
   */

  int Nl = region->level > 0 ? 1 << (region->level - 1) : N/Dimensions.x;
  if (region->n.x == 1 && region->n.y == 1) {
    assert (!GPUContext.fragment_shader);
    csOrigin.x = (region->p.x - X0)/L0*Nl*Dimensions.x;
    csOrigin.y = (region->p.y - Y0)/L0*Nl*Dimensions.x;
    CL_CHECK (clSetKernelArg (shader->kernel, 2, sizeof (struct { int x, y; }), &csOrigin));
    size_t global_work_size[2] = {1, 1};
    size_t local_work_size[2] = {1, 1};
    CL_CHECK (clEnqueueNDRangeKernel (queue, shader->kernel, 2, NULL,
                                      global_work_size, local_work_size,
                                      0, NULL, NULL));
  }

  /**
   * This is a region
   */

  else if (region->n.x || region->n.y) {
    // OpenCL doesn't use glDrawArrays, so we need to use compute shaders
    // This would require a different approach for region rendering
    assert (false);
  }

  else {
    assert (!GPUContext.fragment_shader);
    size_t global_work_size[2] = {
      shader->ng[0]*shader->nwg[0],
      shader->ng[1]*shader->nwg[1]
    };
    size_t local_work_size[2] = {
      shader->nwg[0],
      shader->nwg[1]
    };
    CL_CHECK (clEnqueueNDRangeKernel (queue, shader->kernel, 2, NULL,
                                      global_work_size, local_work_size,
                                      0, NULL, NULL));
  }
  CL_CHECK (clFinish (queue));
  return Nl;
}

void gpu_free_solver (void)
{
  CL_CHECK (clFinish (queue));
  if (queue) CL_CHECK (clReleaseCommandQueue (queue));
  if (ctx) CL_CHECK (clReleaseContext (ctx));
  ctx = NULL;
  queue = NULL;
}

void gpu_synchronize (void)
{
  if (ctx)
    CL_CHECK (clFinish (queue));
}

/* -------------------------------------------------------------------------
   OpenCL Kernel Source String
   Modified to accept an 'offset' parameter and apply it to the input indices.
   ------------------------------------------------------------------------- */
static char kernel_source[] =
"#define REDUCE(reduced,rhs)                                                     \n"
"__kernel void reduce (__global const float* input, __global float* output, const int n, const int offset){\n"
"    __local float sdata[256];                                                   \n"
"    unsigned int tid = get_local_id(0);                                         \n"
"    unsigned int i = get_group_id(0) * get_local_size(0) * 2 + tid;             \n"
"    float reduced =                                                             \n" /* Patched dynamically */
"    if (i < n)                                                                  \n"
"        REDUCE (reduced, input[i + offset]);                                    \n" /* Offset applied here */
"    if (i + get_local_size(0) < n)                                              \n"
"        REDUCE (reduced, input[i + get_local_size(0) + offset]);                \n" /* Offset applied here */
"    sdata[tid] = reduced;                                                       \n"
"    barrier(CLK_LOCAL_MEM_FENCE);                                               \n"
"                                                                                \n"
#if 1 // Completely portable loop: does not depend on a warp size of 32
"    for (unsigned int s = get_local_size(0) / 2; s > 0; s >>= 1) {              \n"
"      if (tid < s) {                                                            \n"
"        REDUCE(sdata[tid], sdata[tid + s]);                                     \n"
"      }                                                                         \n"
"      barrier(CLK_LOCAL_MEM_FENCE); // Forces all hardware types to sync safely \n"
"    }                                                                           \n"
#else // optimised for NVIDIA's warp size of 32
"    for (unsigned int s = get_local_size(0)/2; s > 32; s >>= 1){\n"
"        if (tid < s)                                                            \n"
"            REDUCE (sdata[tid], sdata[tid + s]);                                 \n"
"        barrier(CLK_LOCAL_MEM_FENCE);                                           \n"
"    }                                                                           \n"
"    if (tid < 32) {                                                             \n"
"        __local volatile float* smem = sdata;                                   \n"
"        REDUCE (smem[tid], smem[tid + 32]);                                     \n"
"        REDUCE (smem[tid], smem[tid + 16]);                                     \n"
"        REDUCE (smem[tid], smem[tid + 8]);                                      \n"
"        REDUCE (smem[tid], smem[tid + 4]);                                      \n"
"        REDUCE (smem[tid], smem[tid + 2]);                                      \n"
"        REDUCE (smem[tid], smem[tid + 1]);                                      \n"
"    }                                                                           \n"
#endif
"    if (tid == 0)                                                               \n"
"        output[get_group_id(0)] = sdata[0];                                     \n"
"}                                                                               \n";

/* -------------------------------------------------------------------------
   Dynamic Runtime Kernel Compilation
   ------------------------------------------------------------------------- */
static cl_kernel compile_kernel (const char * start, const char * op)
{
    // 1. Patch the macro definition string
    char * s = kernel_source + strlen("#define REDUCE(reduced,rhs) ");
    memcpy (s, op, strlen (op));
    s += strlen(op); while (*s != '\n') *s++ = ' ';
    
    // 2. Patch the variable initialization line
    s = strstr (kernel_source, "float reduced = ");
    s += strlen ("float reduced = ");
    memcpy (s, start, strlen (start));
    s += strlen(start); while (*s != '\n') *s++ = ' '; 

    // 3. Create OpenCL Program object from the source string
    cl_int err;
    const char* src_ptr = kernel_source;
    size_t src_len = strlen(kernel_source);
    cl_program program = clCreateProgramWithSource(ctx, 1, &src_ptr, &src_len, &err);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "Error creating program from source: %d\n", err);
        exit(EXIT_FAILURE);
    }
    
    // 4. Compile the OpenCL source code
    err = clBuildProgram(program, 0, NULL, "-cl-fast-relaxed-math", NULL, NULL);
    if (err != CL_SUCCESS) {
        size_t log_size;
        cl_device_id device;
        clGetCommandQueueInfo(queue, CL_QUEUE_DEVICE, sizeof(cl_device_id), &device, NULL);
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
        
        char *log = (char *) malloc(log_size);
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, log_size, log, NULL);
        fprintf(stderr, "OpenCL Compilation Failure:\n%s\n", log);
        free(log);
        exit(EXIT_FAILURE);
    }

    // 5. Extract the compiled executable entry point
    cl_kernel kernel = clCreateKernel(program, "reduce", &err);
    if (err != CL_SUCCESS) {
        fprintf(stderr, "Error creating kernel object: %d\n", err);
        exit(EXIT_FAILURE);
    }
    
    clReleaseProgram(program); 
    return kernel;
}

/* -------------------------------------------------------------------------
   Multi-Pass Reduction Execution
   ------------------------------------------------------------------------- */
static float opencl_reduce (cl_mem d_input, const size_t N, const char op, const int initial_offset)
{
    cl_kernel kernel;
    switch (op) {
    case '+': {
        static cl_kernel k = 0;
        if (!k) k = compile_kernel ("0.0f;", "reduced += rhs");
        kernel = k;
        break;
    }
    case 'm': {
        static cl_kernel k = 0;
        if (!k) k = compile_kernel ("3e38f;", "reduced = min(reduced,rhs)");
        kernel = k;
        break;
    }
    case 'M': {
        static cl_kernel k = 0;
        if (!k) k = compile_kernel ("-3e38f;", "reduced = max(reduced,rhs)");
        kernel = k;
        break;
    }
    default:
        assert (0);
    }

    /* Track state allocations globally across sequential steps */
    static size_t Np = 0;
    static cl_mem d_output_a = NULL, d_output_b = NULL;
    cl_int err;

    if (N > Np) {
        const size_t max_blocks = (N + 511)/512;
        if (d_output_a) {
            clReleaseMemObject(d_output_a);
            clReleaseMemObject(d_output_b);
        }
        d_output_a = clCreateBuffer(ctx, CL_MEM_READ_WRITE, max_blocks * sizeof(float), NULL, &err);
        d_output_b = clCreateBuffer(ctx, CL_MEM_READ_WRITE, max_blocks * sizeof(float), NULL, &err);
        Np = N;
    }
        
    /* Execute down-scaling loops */
    cl_mem input = d_input;
    cl_mem output = d_output_a;
    size_t current_n = N;
    
    // The initial pass applies your calculated data offset.
    // Subsequent multi-pass reductions read from intermediate arrays starting at 0.
    int current_offset = initial_offset; 

    while (current_n > 1) {
        const size_t local_size = 256;
        size_t blocks = (current_n + local_size*2 - 1)/(local_size*2);
        size_t global_size = blocks * local_size; 
        int n_param = (int)current_n;

        // Pass arguments sequentially: input, output, n, offset
        clSetKernelArg(kernel, 0, sizeof(cl_mem), &input);
        clSetKernelArg(kernel, 1, sizeof(cl_mem), &output);
        clSetKernelArg(kernel, 2, sizeof(int), &n_param);
        clSetKernelArg(kernel, 3, sizeof(int), &current_offset);

        err = clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &global_size, &local_size, 0, NULL, NULL);
        if (err != CL_SUCCESS) {
            fprintf(stderr, "Error enqueuing reduction kernel stage: %d\n", err);
            exit(EXIT_FAILURE);
        }

        current_n = blocks;
        input = output;
        output = (output == d_output_a) ? d_output_b : d_output_a;
        
        // After the first loop pass, the output is packed cleanly starting at index 0
        current_offset = 0; 
    }

    /* Pull the single final result back to host CPU */
    float gpu_reduce;
    clEnqueueReadBuffer(queue, input, CL_TRUE, 0, sizeof(float), &gpu_reduce, 0, NULL, NULL);

    return gpu_reduce;
}

/* -------------------------------------------------------------------------
   Main Top-Level API Interface Wrapper
   ------------------------------------------------------------------------- */
double gpu_reduction (size_t offset,
                      const char op,
                      const RegionParameters * region,
                      GPUData * data,
                      size_t nb)
{
  if (region->n.x == 1 && region->n.y == 1) {
    int i = (region->p.x - X0)/L0*N;
    int j = (region->p.y - Y0)/L0*N;
    if (i < 0 || i >= N || j < 0 || j >= N)
      return 0.;
    offset += i*N + j;
    nb = 1;
  }
    
  float result;
  if (nb == 1)
    CL_CHECK (clEnqueueReadBuffer (queue, ssbo,
                                   CL_TRUE, // Blocking read
                                   offset*sizeof(float), sizeof(float), &result, 0, NULL, NULL));
  else
    opencl_reduce (ssbo, nb, op, (int)offset); // fixme: what about very large offsets??
  return result;
}
