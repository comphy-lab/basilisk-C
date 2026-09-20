#define dimension 3
#define GRIDNAME "Multigrid 3D (gpu)"
#include "../gpu-multigrid.h"
#pragma autolink -L$BASILISK/grid/gpu -lgpu -lerrors -lglfw -ldl

static void gpu_multigrid3D_methods()
{
  multigrid_methods();
  boundary_level = gpu_boundary_level;
}
