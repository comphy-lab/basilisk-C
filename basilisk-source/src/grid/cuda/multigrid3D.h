#define dimension 3
#define GRIDNAME "Multigrid 3D (cuda)"
#define _CUDA 1
#include "../gpu-multigrid.h"
#pragma autolink -L$BASILISK/grid/cuda -lbuda -lcuda -lnvrtc -L$BASILISK/grid/gpu -lerrors

static void cuda_multigrid3D_methods()
{
  multigrid_methods();
  boundary_level = gpu_boundary_level;
}
