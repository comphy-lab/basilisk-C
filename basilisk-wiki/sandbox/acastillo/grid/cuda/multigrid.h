#define SINGLE_PRECISION 1
#define GRIDNAME "Multigrid (cuda)"
#define _CUDA 1
#include "../gpu-multigrid.h"
#if DOUBLE_PRECISION
# pragma autolink -L$BASILISK_INCLUDE_PATH/acastillo/grid/cuda -lbuda64 -lcuda -lnvrtc -L$BASILISK/grid/gpu -lerrors
#else
# pragma autolink -L$BASILISK_INCLUDE_PATH/acastillo/grid/cuda -lbuda -lcuda -lnvrtc -L$BASILISK/grid/gpu -lerrors
#endif

static void cuda_multigrid_methods()
{
  multigrid_methods();
  boundary_level = gpu_boundary_level;
}
