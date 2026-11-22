#define dimension 3
#define DOUBLE_PRECISION 1
#define GRIDNAME "Multigrid 3D (gpu FP64)"
#include "../gpu-multigrid.h"
#pragma autolink -L$BASILISK_INCLUDE_PATH/acastillo/grid/gpu -lgpu64 -L$BASILISK/grid/gpu -lerrors -lglfw -ldl

static void gpu_multigrid3D_fp64_methods()
{
  multigrid_methods();
  boundary_level = gpu_boundary_level;
}
