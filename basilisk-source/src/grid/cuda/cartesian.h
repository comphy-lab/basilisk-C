#define GRIDNAME "Cartesian (cuda)"
#define _CUDA 1
#include "../gpu-cartesian.h"
#pragma autolink -L$BASILISK/grid/cuda -lbuda -lcuda -lnvrtc -L$BASILISK/grid/gpu -lerrors
               
static void cuda_cartesian_methods()
{
  cartesian_methods();
  boundary_level = gpu_boundary_level;
}
