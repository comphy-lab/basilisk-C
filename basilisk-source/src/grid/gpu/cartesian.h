#define GRIDNAME "Cartesian (GPU)"
#include "../gpu-cartesian.h"
#pragma autolink -L$BASILISK/grid/gpu -lgpu -lerrors -lglfw -ldl

static void gpu_cartesian_methods()
{
  cartesian_methods();
  boundary_level = gpu_boundary_level;
}
