/**
# Evaporation

The evaporation model is a diffusion-driven multicomponent phase change model,
which can resolve multiple chemical species in liquid and gas phases, immersed
in a non-isothermal environment. The current module simplifies setting up
evaporation cases by automatically adjusting the `pcm` (phase change model)
parameters.
*/

#include "phasechange.h"
#include "multicomponent.h"

// fixme: check compatibility with restore
event init (i = 0) {
  pcm.shifting = SHIFT_TO_LIQUID;
  pcm.diffusion = EXPLICIT_IMPLICIT;
  if (nv == 1)
    pcm.byrhogas = true;
}
