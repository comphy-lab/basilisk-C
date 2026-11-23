#ifndef SINGLE_PRECISION
# ifdef DOUBLE_PRECISION
#  define SINGLE_PRECISION (!DOUBLE_PRECISION)
# else
#  define SINGLE_PRECISION 0
# endif
#endif
#ifndef DOUBLE_PRECISION
# define DOUBLE_PRECISION (!SINGLE_PRECISION)
#endif
#define _GPU 1
#define GRIDPARENT Multigrid
#define shift_level(d) (multigrid->shift[d])
#define field_size() (multigrid->shift[depth() + 1])
#define grid_data() (multigrid->d)
#define field_offset(s, level) (shift_level((level) ? (level) - 1 : depth()) + (s).i*field_size())

#if _CUDA
void cuda_set_device (int local_rank);
void cuda_halo_sync (void);
int cuda_query_mpi_support (void);
size_t cuda_halo_pack (const int * fields, int nfields, size_t fsize, size_t shift,
		       int ny_tot, int nz_tot, int i0, int j0, int k0,
		       int ci, int cj, int ck, int slot);
void * cuda_halo_snd_ptr (int slot, int to_host, size_t nbytes);
void * cuda_halo_rcv_ptr (int slot, int from_host, size_t nbytes);
void cuda_halo_unpack (const int * fields, int nfields, size_t fsize, size_t shift,
		       int ny_tot, int nz_tot, int i0, int j0, int k0,
		       int ci, int cj, int ck, int slot, int from_host, size_t nbytes);
#endif

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
#include "grid/stencils.h"
#include "gpu/gpu.h"
#include "grid/multigrid-common.h"
#include "grid/gpu/backend.h"
#if _MPI
# include "multigrid-mpi.h"
#endif

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

static inline void gpu_subtree_size (scalar size, bool leaves)
{
  subtree_size (size, leaves);
  if (size.gpu.stored < 0) {
    gpu_cpu_sync_scalar (size.i, size.block, grid_data(), field_size(), GPU_READ);
    size.gpu.stored = 0;
  }
}
#define subtree_size(s, l) gpu_subtree_size(s, l)
