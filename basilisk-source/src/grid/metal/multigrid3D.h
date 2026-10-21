#define dimension 3
#define GRIDNAME "Multigrid 3D (metal)"
#define _CUDA 1
#define _METAL 1
#include "../gpu-multigrid.h"
#pragma autolink -L$BASILISK/grid/metal -lmetal -framework Metal -framework Foundation -L$BASILISK/grid/gpu -lerrors

/**
The Metal Shading Language does not support double-precision floating
point numbers. */

#if !SINGLE_PRECISION
# error Metal does not support double precision
#endif

static void metal_multigrid3D_methods()
{
  multigrid_methods();
  boundary_level = gpu_boundary_level;
}
