#define dimension 3
#define SINGLE_PRECISION 1
#define GRIDNAME "Multigrid 3D (gpu)"
#include "../gpu-multigrid.h"
#if DOUBLE_PRECISION
# pragma autolink -L$BASILISK_INCLUDE_PATH/acastillo/grid/gpu -lgpu64 -L$BASILISK/grid/gpu -lerrors -lglfw -ldl
#else
# pragma autolink -L$BASILISK_INCLUDE_PATH/acastillo/grid/gpu -lgpu -L$BASILISK/grid/gpu -lerrors -lglfw -ldl
#endif

static void gpu_multigrid3D_methods()
{
  multigrid_methods();
  boundary_level = gpu_boundary_level;
}
