#define GRIDNAME "Multigrid (opencl)"
#define _CUDA 1
#define _OPENCL 1
#include "../gpu-multigrid.h"
#pragma autolink -L$BASILISK/grid/opencl -locl -lOpenCL -L$BASILISK/grid/gpu -lerrors

static void opencl_multigrid_methods()
{
  multigrid_methods();
  boundary_level = gpu_boundary_level;
}
