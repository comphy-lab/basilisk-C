#define GRIDNAME "Cartesian (hip)"
#define _CUDA 1
#include "../gpu-cartesian.h"

/**
By default _HIP_AMD is not defined and the CUDA version is used. To
use the AMD version compile your Basilisk program with -D_HIP_AMD=1. */

#if _HIP_AMD
  #pragma autolink -L$BASILISK/grid/hip -lhipamd -lamdhip64 -lhiprtc -L$BASILISK/grid/gpu -lerrors
#else // HIP-on-NVIDIA
  #pragma autolink -L$BASILISK/grid/hip -lhipcuda -lcuda -lcudart -lnvrtc -L$BASILISK/grid/gpu -lerrors
#endif

static void hip_cartesian_methods()
{
  cartesian_methods();
  boundary_level = gpu_boundary_level;
}
