#define DOUBLE_PRECISION 1
#define GRIDNAME "Multigrid (gpu FP64)"
#include "../gpu-multigrid.h"
#pragma autolink -L$BASILISK_INCLUDE_PATH/acastillo/grid/gpu -lgpu64 -L$BASILISK/grid/gpu -lerrors -lglfw -ldl

static void gpu_multigrid_fp64_methods()
{
  multigrid_methods();
  boundary_level = gpu_boundary_level;
}
