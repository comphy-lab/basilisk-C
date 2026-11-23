#define dimension 3
#define DOUBLE_PRECISION 1
#define GRIDNAME "Multigrid 3D (cuda FP64)"
#define _CUDA 1
#include "../gpu-multigrid.h"
#pragma autolink -L$BASILISK_INCLUDE_PATH/acastillo/grid/cuda -lbuda64 -lcuda -lnvrtc -ldl -L$BASILISK/grid/gpu -lerrors

static void cuda_multigrid3D_fp64_methods()
{
  multigrid_methods();
  boundary_level = gpu_boundary_level;
}
