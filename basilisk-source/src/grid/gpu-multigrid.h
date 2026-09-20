#if DOUBLE_PRECISION
# define SINGLE_PRECISION 0
#else
# define SINGLE_PRECISION 1
#endif
#define _GPU 1
#define GRIDPARENT Multigrid
#define shift_level(d) (multigrid->shift[d])
#define field_size() (multigrid->shift[depth() + 1])
#define grid_data() (multigrid->d)
#define field_offset(s, level) (shift_level((level) ? (level) - 1 : depth()) + (s).i*field_size())

#if dimension == 2
#define GPU_CODE()							\
  "#define valt(s,di,dj,dk)"						\
  "  _data_val(_index(s,dk), point.j + (dj) + "				\
  " (point.i + (di))*(point.n.y + 2*GHOSTS) + _shift[point.level])\n"	\
  "#define val_red_(s) _data_val((s).i, point.j - GHOSTS +"		\
  "  (point.i - GHOSTS)*NY + _shift[point.level])\n"			\
  "#define fine(a,di,dj,dk)"						\
  "  _data_val(_index(a,dk), 2*point.j - GHOSTS + (dj) +"               \
  "        (2*point.i - GHOSTS + (di))*(point.n.y*2 + 2*GHOSTS) +"	\
  "        _shift[point.level + 1])\n"					\
  "#define coarse(a,di,dj,dk)"						\
  "  _data_val(_index(a,dk), (point.j + GHOSTS)/2 + (dj) +"		\
  "        ((point.i + GHOSTS)/2 + (di))*(point.n.y/2 + 2*GHOSTS) +"	\
  "        _shift[point.level - 1])\n"
#elif dimension == 3
#define GPU_CODE()							\
  "#define valt(s,di,dj,dk)"						\
  "  _data_val(_index(s,0), point.k + (dk) +"				\
  " (point.n.z + 2*GHOSTS)*(point.j + (dj) +"                           \
  " (point.i + (di))*(point.n.y + 2*GHOSTS)) +"                         \
  " _shift[point.level])\n"						\
  "#define val_red_(s) _data_val((s).i, point.k - GHOSTS +"		\
  " NZ*(point.j - GHOSTS + (point.i - GHOSTS)*NY) +"			\
  " _shift[point.level])\n"						\
  "#define fine(a,di,dj,dk)"						\
  "  _data_val(_index(a,0), 2*point.k - GHOSTS + (dk) +"		\
  " (point.n.z*2 + 2*GHOSTS)*(2*point.j - GHOSTS + (dj) +"		\
  " (2*point.i - GHOSTS + (di))*(point.n.y*2 + 2*GHOSTS)) +"            \
  " _shift[point.level + 1])\n"                                         \
  "#define coarse(a,di,dj,dk)"						\
  "  _data_val(_index(a,0), (point.k + GHOSTS)/2 + (dk) +"		\
  " (point.n.z/2 + 2*GHOSTS)*((point.j + GHOSTS)/2 + (dj) +"            \
  " ((point.i + GHOSTS)/2 + (di))*(point.n.y/2 + 2*GHOSTS)) +"          \
  " _shift[point.level - 1])\n"
#endif

static bool _gpu_done_ = false;

#include "multigrid.h"
#include "stencils.h"
#include "gpu/gpu.h"
#include "multigrid-common.h"
#include "gpu/backend.h"

void realloc_scalar_gpu (int size)
{
  for (scalar s in baseblock) {
    if (s.gpu.stored < 0)
      gpu_cpu_sync_scalar (s.i, s.block, grid_data(), field_size(), GPU_READ);
    s.gpu.stored = 1; // only stored on the CPU
  }
  realloc_scalar_cpu (size);
  realloc_ssbo (field_size());
}

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
