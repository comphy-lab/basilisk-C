#define GRIDNAME "Multigrid (hip)"
#define _CUDA 1
#include "../gpu-multigrid.h"

/**
Note that the libraries `-lcuda -lcudart -lnvrtc` are only required for the
HIP-on-CUDA version. This will need to be adapted for the HIP/ROCm
version. See also the [Makefile]() for compilation of the
corresponding `libhip.a` library. */

#pragma autolink -L$BASILISK/grid/hip -lhip -lcuda -lcudart -lnvrtc -L$BASILISK/grid/gpu -lerrors

static void hip_multigrid_methods()
{
  multigrid_methods();
  boundary_level = gpu_boundary_level;
}
