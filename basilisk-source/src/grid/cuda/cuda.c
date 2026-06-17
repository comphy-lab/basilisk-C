/**
# CUDA backend for GPUs

This relies on the CUDA and
[NVRTC](https://docs.nvidia.com/cuda/nvrtc/index.html) libraries which
can be installed on Debian systems using:

~~~bash
apt install nvidia-driver nvidia-cuda-dev
~~~

For some reason the performances are very low compared to the [GLSL
backend](/src/grid/gpu/opengl.c) (on the same nvidia card). */

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "a32.h"

typedef struct { double x, y, z; } coord;
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

typedef struct _External External;

struct _External {
  char * name;    // the name of the variable
  void * pointer; // a pointer to the data
  int type;       // the type of the variable
  int nd;         // the number of pointer dereferences or attribute offset or enum constant
  char reduct;    // the reduction operation
  char global;    // is it a global variable?
  char constant;  // is it a constant?
  void * data;    // the dimensions (int *) for arrays or the code (char *) for functions
  scalar s;       // used for reductions on GPUs
  External * externals, * next;
  int used;
};

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

#include "../gpu/backend.h"

#include <cuda.h>
#include <nvrtc.h>

static CUdeviceptr ssbo = 0;
static CUdevice dev = 0;
static CUcontext ctx = 0;
static CUstream stream = 0;

#define CUDA_CHECK(x)                                                   \
  do {                                                                  \
    CUresult err = (x);                                                 \
    if (err != CUDA_SUCCESS) {                                          \
      const char *msg;                                                  \
      cuGetErrorName(err, &msg);                                        \
      fprintf(stderr, "%s:%d: CUDA error: %s\n", __FILE__, __LINE__, msg); \
      exit(1);                                                          \
    }                                                                   \
  } while(0)

#define NVRTC_CHECK(x)                          \
  do {                                          \
    nvrtcResult err = (x);                      \
    if (err != NVRTC_SUCCESS) {                 \
      fprintf(stderr, "%s:%d: NVRTC error: %s\n", __FILE__, __LINE__,   \
              nvrtcGetErrorString(err));        \
      exit(1);                                  \
    }                                           \
  } while(0)

typedef struct {
  CUdeviceptr location;
  int type, nd, local;
  size_t size;
  void * pointer;
} MyUniform;

struct _Shader {
  unsigned ng[2], nwg[2];
  CUdeviceptr _data, csOrigin;
  MyUniform * uniforms;
  CUmodule module;
  CUfunction kernel;
};

void free_shader (Shader * s)
{
  free (s->uniforms);
  free (s);
}

static
void architecture (char * arch)
{
  int major, minor;
  CUdevice cuDevice;
  CUDA_CHECK (cuDeviceGet(&cuDevice, 0));
  CUDA_CHECK (cuDeviceGetAttribute (&major, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR,
                                    cuDevice));
  CUDA_CHECK (cuDeviceGetAttribute (&minor, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR,
                                    cuDevice));
  sprintf (arch, "--gpu-architecture=sm_%d%d", major, minor);
  // fixme: not sure whether this should be compute_%d%d or sm_%d%d ??
}

Shader * load_normal_shader (const char * fs, const char * func, const char * file, int line)
{
  //  fputs (fs, stderr);
  nvrtcProgram prog;
  NVRTC_CHECK(nvrtcCreateProgram (&prog, fs,
                                  "kernel.cu",
                                  0,
                                  NULL,
                                  NULL
                                  ));
  char arch[] = "--gpu-architecture=compute_????";
  architecture (arch);
  const char *opts[] = {
    "--std=c++11",
    arch,
    "-default-device",
    "-diag-suppress=177",
    "--ptxas-options=-O3",
    "--extra-device-vectorization",
    "--restrict",
    "-use_fast_math",
  };

  nvrtcResult compile_res = nvrtcCompileProgram (prog, sizeof(opts)/sizeof(char *), opts);

  if (compile_res != NVRTC_SUCCESS) {
    size_t logSize;
    NVRTC_CHECK (nvrtcGetProgramLogSize (prog, &logSize));
    if (logSize > 1) {
      char * log = (char *) malloc (logSize);
      NVRTC_CHECK (nvrtcGetProgramLog (prog, log));
      // fputs (log, stderr);
      char * error = gpu_errors (log, fs, NULL, "CUDA");
      fputs (error, stderr);
      free (error);
      free (log);
    }
    return NULL;
  }

  // ------------------------------------------------------------
  // Get PTX
  // ------------------------------------------------------------

  size_t ptxSize;
  NVRTC_CHECK (nvrtcGetPTXSize (prog, &ptxSize));
  char ptx[ptxSize];
  NVRTC_CHECK (nvrtcGetPTX (prog, ptx));
  NVRTC_CHECK (nvrtcDestroyProgram (&prog));

#if SINGLE_PRECISION
  //  fputs (ptx, stderr);
  if (strstr (ptx, ".f64"))
    fprintf (stderr, "%s:%d: warning: CUDA: found FP64 assembly in single precision mode\n",
             file, line);
#endif // SINGLE_PRECISION
  
  // ------------------------------------------------------------
  // Load PTX module
  // ------------------------------------------------------------

  Shader * shader = calloc (1, sizeof (Shader));
  CUDA_CHECK (cuModuleLoadData (&shader->module, ptx));
  CUDA_CHECK (cuModuleGetFunction (&shader->kernel, shader->module, func));
  return shader;
}

bool gpu_init_context (GPUData ** data)
{
  bool initialized = ctx;
  if (!initialized) {
    CUDA_CHECK (cuInit (0));
    CUDA_CHECK (cuDeviceGet (&dev, 0));
    CUDA_CHECK (cuCtxCreate (&ctx, 0, dev));
  }
  *data = NULL;
  return !initialized;
}

void gpu_free_context (GPUData * data)
{
  if (ssbo) {
    CUDA_CHECK (cuMemFree (ssbo));
    ssbo = 0;
  }
  GPUContext.current_size = 0;
}

void realloc_ssbo (size_t field_size)
{
  if (!datasize)
    return;
  size_t totalsize = field_size*datasize;
  assert (totalsize > GPUContext.current_size);
  CUdeviceptr ptr;
  CUDA_CHECK (cuMemAlloc (&ptr, totalsize)); // fixme: allocates memory twice
  if (GPUContext.current_size > 0) {
    CUDA_CHECK (cuMemcpyDtoD (ptr, ssbo, GPUContext.current_size));
    CUDA_CHECK (cuMemFree (ssbo));
  }
  ssbo = ptr;
  GPUContext.current_size = totalsize;
}

void gpu_cpu_sync_scalar (int i, int block, char * data, size_t field_size, SyncMode mode)
{
  size_t size = field_size*sizeof(real), offset = i*size, totalsize = block*size;
  char * cd = data + offset;
  CUdeviceptr gd = ssbo + offset;
  if (mode == GPU_READ)
    CUDA_CHECK (cuMemcpyDtoH (cd, gd, totalsize));
  else if (mode == GPU_WRITE)
    CUDA_CHECK (cuMemcpyHtoD (gd, cd, totalsize));
  else
    assert (false);
}

void reset_scalar (int i, int block, size_t field_size, double val)
{
  size_t size = field_size*sizeof(real);
  size_t offset = i*size, totalsize = max(block, 1)*size;
#if SINGLE_PRECISION
  float fval = val;
  uint32_t bits;
  memcpy (&bits, &fval, sizeof(bits));
  CUDA_CHECK (cuMemsetD32 (ssbo + offset, bits, totalsize/sizeof(float)));
#else
  fprintf (stderr, "%s:%d: error: not implemented yet\n");
#endif
}

void finalize_shader (Shader * shader, External * externals, External * merged,
                      unsigned ng[2], unsigned nwg[2])
{
  for (int i = 0; i < 2; i++)
    shader->ng[i] = ng[i], shader->nwg[i] = nwg[i];
  
  /**
  Get the SSBO pointer. */

  size_t size;
  CUDA_CHECK (cuModuleGetGlobal(&shader->_data, &size, shader->module, "_data"));
  assert (size == sizeof (ssbo));

  /**
  Get the csOrigin pointer. */

  CUDA_CHECK (cuModuleGetGlobal(&shader->csOrigin, &size, shader->module, "csOrigin"));
  assert (size == 2*sizeof (int));
  
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
    if (IS_EXTERNAL_CONSTANT(g)) continue;
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
      CUdeviceptr location = 0;
      size_t size;
      cuModuleGetGlobal (&location, &size, shader->module, name);
      if (location) {
        //        fprintf (stderr, "%s %d %ld\n", name, g->type, size);
	// not an array or just a one-dimensional array
	assert (!g->nd);
	assert (!g->data || ((int *)g->data)[1] == 0);
	int nd = g->data ? ((int *)g->data)[0] : 1;
        switch (g->type) {
        case sym_INT: case sym_LONG:
          assert (size == sizeof(int)*nd); break;
        case sym_FLOAT:
          assert (size == sizeof(float)*nd); break;
        case sym_VEC4:
          nd *= 4;
          assert (size == sizeof(float)*nd); break;
        case sym_BOOL:
          assert (size == sizeof(bool)*nd); break;
#if SINGLE_PRECISION
        case sym_DOUBLE:
          assert (size == sizeof(float)*nd); break;
        case sym__COORD:
          nd *= 2;
          assert (size == sizeof(float)*nd); break;
        case sym_COORD:
          nd *= 3;
          assert (size == sizeof(float)*nd); break;
#else // DOUBLE_PRECISION
        case sym_DOUBLE:
          assert (size == sizeof(double)*nd); break;
        case sym__COORD:
          nd *= 2;
          assert (size == sizeof(double)*nd); break;
        case sym_COORD:
          nd *= 3;
          assert (size == sizeof(double)*nd); break;
#endif // DOUBLE_PRECISION
        default:
          assert (false);
        }
	shader->uniforms = realloc (shader->uniforms, (nuniforms + 2)*sizeof(MyUniform));
	shader->uniforms[nuniforms] = (MyUniform){
	  .location = location, .type = g->type, .nd = nd, .size = size,
	  .local = g->global == 1 ? -1 : g->used - 1,
	  .pointer = g->global == 1 ? g->pointer : NULL };
	shader->uniforms[nuniforms + 1].type = 0;
	nuniforms++;
	// uniforms refering to local variables must be in the 'externals' local list
	assert (g->global == 1 || g->used);
      }
      else
        fprintf (stderr, "%s not found\n", name);
      free (name);
    }
  }
}

void post_setup_shader (Shader * shader, External * externals)
{
  
  /**
  Set SSBO pointer. */
  
  assert (ssbo);
  CUDA_CHECK (cuMemcpyHtoDAsync (shader->_data, &ssbo, sizeof (ssbo), stream));

  /**
  ## Set uniforms */

  for (const MyUniform * g = shader->uniforms; g && g->type; g++) {
    void * pointer = g->pointer;
    if (!pointer) {
      assert (g->local >= 0);
      pointer = externals[g->local].pointer;
    }
    switch (g->type) {
    case sym_INT: case sym_FLOAT: case sym_VEC4: case sym_BOOL:
      CUDA_CHECK (cuMemcpyHtoDAsync (g->location, pointer, g->size, stream));
      break;
    case sym_LONG: {
      int p[g->nd];
      long * data = pointer;
      for (int i = 0; i < g->nd; i++)
	p[i] = data[i];
      CUDA_CHECK (cuMemcpyHtoDAsync (g->location, p, g->size, stream));
      break;
    }
#if SINGLE_PRECISION
    case sym_DOUBLE: case sym__COORD: case sym_COORD: {
      float p[g->nd];
      double * data = pointer;
      for (int i = 0; i < g->nd; i++)
	p[i] = data[i];
      CUDA_CHECK (cuMemcpyHtoDAsync (g->location, p, g->size, stream));
      break;
    }
#else // DOUBLE_PRECISION
    case sym_DOUBLE: case sym__COORD: case sym_COORD:
      CUDA_CHECK (cuMemcpyHtoDAsync (g->location, pointer, g->size, stream));
      break;
#endif // DOUBLE_PRECISION
    default:
      assert (false);
    }
  }
}

int run_shader (const Shader * shader, const RegionParameters * region)
{
  /**
  ## Render 

  If this is a `foreach_point()` iteration, we access a single point */

  int Nl = region->level > 0 ? 1 << (region->level - 1) : N/Dimensions.x;
  if (region->n.x == 1 && region->n.y == 1) {
    int csOrigin[] = {
      (region->p.x - X0)/L0*Nl*Dimensions.x,
      (region->p.y - Y0)/L0*Nl*Dimensions.x
    };
    CUDA_CHECK (cuMemcpyHtoDAsync (shader->csOrigin, csOrigin, 2*sizeof(int), stream));
    assert (!GPUContext.fragment_shader);
    CUDA_CHECK (cuLaunchKernel (shader->kernel,
                                1, 1, 1,
                                1, 1, 1,
                                0, stream, NULL, NULL));
  }

  /**
  This is a region */
  
  else if (region->n.x || region->n.y) {
#if 1
    assert (false);
#else
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
#endif
  }

  else {
    assert (!GPUContext.fragment_shader);
    CUDA_CHECK (cuLaunchKernel (shader->kernel,
                                shader->ng[0], shader->ng[1], 1,
                                shader->nwg[0], shader->nwg[1], 1,
                                0, stream, NULL, NULL));
  }
  return Nl;
}

void gpu_free_solver (void)
{
  CUDA_CHECK (cuCtxSynchronize ());
  CUDA_CHECK (cuCtxDestroy (ctx));
  ctx = NULL;
}

void gpu_synchronize()
{
  if (ctx)
    CUDA_CHECK (cuCtxSynchronize ());
}

/**
## Reductions */
   
static char kernel_source[] =
"#define REDUCE(reduced,rhs) reduced += rhs                                          \n"
"extern \"C\"\n"
"__global__ void reduce (const float* input, float* output, int n){\n"
"    __shared__ float sdata[256];                       \n"
"    unsigned int tid = threadIdx.x;                    \n"
"    unsigned int i = blockIdx.x * blockDim.x * 2 + tid;\n"
"    float reduced = 0.0f;                              \n"
"    if (i < n)                                         \n"
"        REDUCE (reduced, input[i]);                    \n"
"    if (i + blockDim.x < n)                            \n"
"        REDUCE (reduced, input[i + blockDim.x]);       \n"
"    sdata[tid] = reduced;                              \n"
"    __syncthreads();                                   \n"
"                                                       \n"
"    for (unsigned int s = blockDim.x/2; s > 32; s >>= 1){\n"
"        if (tid < s)                                    \n"
"            REDUCE (sdata[tid], sdata[tid + s]);        \n"
"        __syncthreads();                                \n"
"    }                                                   \n"
"    if (tid < 32) {                                     \n"
"        volatile float* smem = sdata;                   \n"
"        REDUCE (smem[tid], smem[tid + 32]);             \n"
"        REDUCE (smem[tid], smem[tid + 16]);             \n"
"        REDUCE (smem[tid], smem[tid + 8]);              \n"
"        REDUCE (smem[tid], smem[tid + 4]);              \n"
"        REDUCE (smem[tid], smem[tid + 2]);              \n"
"        REDUCE (smem[tid], smem[tid + 1]);              \n"
"    }                                                   \n"
"    if (tid == 0)                                       \n"
"        output[blockIdx.x] = sdata[0];                  \n"
"}                                                       \n";

static CUfunction compile_kernel (const char * start, const char * op)
{
  nvrtcProgram prog;
  char * s = kernel_source + strlen("#define REDUCE(reduced,rhs) ");
  memcpy (s, op, strlen (op));
  s += strlen(op); while (*s != '\n') *s++ = ' ';
  s = strstr (kernel_source, "float reduced = ");
  s += strlen ("float reduced = ");
  memcpy (s, start, strlen (start));
  s += strlen(start); while (*s != '\n') *s++ = ' '; 
  NVRTC_CHECK (nvrtcCreateProgram (&prog, kernel_source, "reduce.cu",
                                   0, NULL, NULL));
  char arch[] = "--gpu-architecture=compute_????";
  architecture (arch);
  const char* options[] = {
    arch,
    "--std=c++11",
    "--use_fast_math"
  };
  if (nvrtcCompileProgram (prog, 3, options) != NVRTC_SUCCESS) {
    fputs (kernel_source, stderr);
    size_t log_size;
    NVRTC_CHECK (nvrtcGetProgramLogSize (prog, &log_size));
    if (log_size) {
      char * log = (char *)malloc (log_size);
      NVRTC_CHECK (nvrtcGetProgramLog (prog, log));
      fprintf (stderr, "%s:%d: %s\n", __FILE__, __LINE__, log);
      free (log);
    }
    exit (EXIT_FAILURE);
  }

  /*    Extract PTX    */
  size_t ptx_size;
  NVRTC_CHECK (nvrtcGetPTXSize (prog, &ptx_size));
  char * ptx = (char *) malloc (ptx_size);
  NVRTC_CHECK (nvrtcGetPTX (prog, ptx));
  NVRTC_CHECK (nvrtcDestroyProgram (&prog));

  /*    Load module    */
  CUmodule module;
  CUDA_CHECK (cuModuleLoadData (&module, ptx));
  free (ptx);

  CUfunction kernel;
  CUDA_CHECK(cuModuleGetFunction (&kernel, module, "reduce"));
  return kernel;
}

static
float cuda_reduce (CUdeviceptr d_input, const size_t N, const char op)
{
  /*    Compile kernel with NVRTC    */
  
  CUfunction kernel;
  switch (op) {
  case '+': {
    static CUfunction k = 0;
    if (!k) k = compile_kernel ("0.0f;", "reduced += rhs");
    kernel = k;
    break;
  }
  case 'm': {
    static CUfunction k = 0;
    if (!k) k = compile_kernel ("3e38f;", "reduced = min(reduced,rhs)");
    kernel = k;
    break;
  }
  case 'M': {
    static CUfunction k = 0;
    if (!k) k = compile_kernel ("-3e38f;", "reduced = max(reduced,rhs)");
    kernel = k;
    break;
  }
  default:
    assert (0);
  }

  /*    Allocate device memory    */

  static size_t Np = 0;
  static CUdeviceptr d_output_a = 0, d_output_b = 0;
  if (N > Np) {
    const size_t max_blocks = (N + 511)/512;
    if (d_output_a) {
      CUDA_CHECK (cuMemFree (d_output_a));
      CUDA_CHECK (cuMemFree (d_output_b));
    }
    CUDA_CHECK (cuMemAlloc (&d_output_a, max_blocks*sizeof(float)));
    CUDA_CHECK (cuMemAlloc (&d_output_b, max_blocks*sizeof(float)));
    Np = N;
  }
    
  /*    Multi-pass reduction    */

  CUdeviceptr input = d_input, output = d_output_a;
  size_t current_n = N;
  while (current_n > 1) {
    const size_t threads = 256;
    size_t blocks = (current_n + threads*2 - 1)/(threads*2);
    void * args[] = {&input, &output, &current_n};
    CUDA_CHECK (cuLaunchKernel (kernel,
                                blocks, 1, 1, threads,
                                1, 1, 0, stream,
                                args, NULL));
    //    CUDA_CHECK (cuCtxSynchronize());
    current_n = blocks;
    input = output;
    output = (output == d_output_a) ? d_output_b : d_output_a;
  }

  /*        Copy result back        */  
  float gpu_reduce;
  CUDA_CHECK (cuMemcpyDtoH (&gpu_reduce, input, sizeof(float)));
#if 0
  /*    Cleanup    */
  CUDA_CHECK(cuMemFree(d_output_a));
  CUDA_CHECK(cuMemFree(d_output_b));
  CUDA_CHECK(cuModuleUnload(module));
#endif
  return gpu_reduce;
}

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
  return cuda_reduce (ssbo + offset*sizeof(real), nb, op);
}
