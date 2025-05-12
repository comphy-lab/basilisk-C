/** 
When we want to use the conserving method solver, on an axisymetric case, we
need to include "axi.h" after "navier-stokes/conserving.h"*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"

int main () {
  init_grid(1);
  run();
}
