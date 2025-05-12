#if DOUBLE_PRECISION
# define SINGLE_PRECISION 0
#else
# define SINGLE_PRECISION 1
#endif
#define _GPU 1
#define GRIDNAME "Multigrid (GPU)"
#define GRIDPARENT Multigrid
#define shift_level(d) (multigrid->shift[d])
#define field_size() (multigrid->shift[depth() + 1])
#define grid_data() (multigrid->d)
#define field_offset(s, level) (shift_level((level) ? (level) - 1 : depth()) + (s).i*field_size())

#define GPU_CODE()							\
  "#define valt(s,k,l,m)"						\
  "  _data_val(_index(s,m), point.j + (l) + "				\
  " (point.i + (k))*(point.n.y + 2*GHOSTS) + _shift[point.level])\n"	\
  "#define val_red_(s) _data_val((s).i, point.j - GHOSTS +"		\
  "  (point.i - GHOSTS)*NY + _shift[point.level])\n"			\
  "#define fine(a,k,l,m)"						\
  "  _data_val(_index(a,m), 2*point.j - GHOSTS + (l) +"			\
  "        (2*point.i - GHOSTS + (k))*(point.n.y*2 + 2*GHOSTS) +"	\
  "        _shift[point.level + 1])\n"					\
  "#define coarse(a,k,l,m)"						\
  "  _data_val(_index(a,m), (point.j + GHOSTS)/2 + (l) +"		\
  "        ((point.i + GHOSTS)/2 + (k))*(point.n.y/2 + 2*GHOSTS) +"	\
  "        _shift[point.level - 1])\n"

static bool _gpu_done_ = false;

#include "../multigrid.h"
#include "../stencils.h"
#include "gpu.h"
#include "../multigrid-common.h"

static void gpu_multigrid_methods()
{
  multigrid_methods();
  boundary_level = gpu_boundary_level;
}

#include "grid.h"
