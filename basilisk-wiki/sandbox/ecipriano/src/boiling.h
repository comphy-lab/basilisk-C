/**
# Boiling

The boiling model is a temperature-gradient-driven phase change model, which
assumes the interface to be at the saturation temperature. This module
simplifies setting up boiling simulations by automatically adding the `pcm`
(phase change model) parameters to the best configuration for boiling.
*/

#include "phasechange.h"
#include "temperature-gradient.h"

// fixme: check compatibility with restore
event init (i = 0) {
  pcm.shifting = SHIFT_TO_GAS;
  pcm.boiling = true;
  pcm.byrhogas = false;
}

