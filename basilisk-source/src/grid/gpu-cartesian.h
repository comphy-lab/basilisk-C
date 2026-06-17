#if DOUBLE_PRECISION
# define SINGLE_PRECISION 0
#else
# define SINGLE_PRECISION 1
#endif
#define _GPU 1
#define GRIDPARENT Cartesian
#define field_size() sq((size_t)N + 2)
#define grid_data() (cartesian->d)
#define field_offset(s,level) ((s).i*field_size())
#define depth() 0
#define GPU_CODE()							\
  "#define valt(s,k,l,m) "						\
  "_data_val(_index(s,m), (point.i + (k))*(N + 2) + point.j + (l))\n" \
  "#define val_red_(s) _data_val((s).i, (point.i - 1)*NY + point.j - 1)\n"

static bool _gpu_done_ = false;

#include "cartesian.h"
@define neighborp(k,l,o) neighbor(k,l,o)
#include "stencils.h"
#include "gpu/gpu.h"
#include "cartesian-common.h"
#include "gpu/backend.h"

typedef struct {
  GRIDPARENT parent;
  khash_t(INT) * shaders;
  GPUData * data;
} GridGPU;

@ifndef tracing
  @ def tracing(func, file, line) do {
    gpu_synchronize();
    tracing(func, file, line);
  } while(0) @
  @ def end_tracing(func, file, line) do {
    gpu_synchronize();
    end_tracing(func, file, line);
  } while(0) @
@endif

#include "gpu/grid.h"
