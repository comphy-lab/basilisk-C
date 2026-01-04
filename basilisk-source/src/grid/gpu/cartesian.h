#if DOUBLE_PRECISION
# define SINGLE_PRECISION 0
#else
# define SINGLE_PRECISION 1
#endif
#define _GPU 1
#define GRIDNAME "Cartesian (GPU)"
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

#include "../cartesian.h"
@define neighborp(k,l,o) neighbor(k,l,o)
#include "../stencils.h"
#include "gpu.h"
#include "../cartesian-common.h"

static void gpu_cartesian_methods()
{
  cartesian_methods();
  boundary_level = gpu_boundary_level;
}

#include "grid.h"
