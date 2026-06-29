#define GRIDNAME "Cartesian (opencl)"
#define _CUDA 1
#define _OPENCL 1
#include "../gpu-cartesian.h"
#pragma autolink -L$BASILISK/grid/opencl -locl -lOpenCL -L$BASILISK/grid/gpu -lerrors
               
static void opencl_cartesian_methods()
{
  cartesian_methods();
  boundary_level = gpu_boundary_level;
}
