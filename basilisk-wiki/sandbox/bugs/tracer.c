/**
This possible bug is realted to the boundary of a tracer in an axisymmetric simulation with mask().
In the init event the command boundary({T}) trigger the following error:
/basilisk/src/grid/multigrid-common.h:39:error: Program received signal SIGFPE, Arithmetic exception.

A possible workaround is to use foreach_boundary() instead of boundary() for the tracer */
#include "axi.h" 
#include "navier-stokes/centered.h"
#include "tracer.h"

#define TRASH 1

#define LEVEL 7
#define MAXTIME 10

scalar T[];
scalar * tracers = {T};

int main() {
  size(1.);
  init_grid(1 << LEVEL);

  TOLERANCE = 1e-6;

  run();
}

event init(t = 0) {
  /**
  Define a rectangular domain using mask() */
  mask (y > 0.5 ? top : none);
  
  /**
  Initialise tracer */
  foreach()
    T[] = x;
  
  /**
  Apply boundary condition (using boundary({T}), you get T = nan at the boundaries)  but
  adding cm to the boundary remedies the bug, */
#if 1  
  boundary({cm, T})
#else
  foreach_boundary(bottom)
    T[ghost] = x;
  foreach_boundary(top)
    T[ghost] = x;
  foreach_boundary(left)
    T[ghost] = x;
  foreach_boundary(right)
    T[ghost] = x;
#endif
}

/**
Simulation snapshot at the end */
event snapshot (t = MAXTIME) {
  output_ppm(T, file = "T.png");
}