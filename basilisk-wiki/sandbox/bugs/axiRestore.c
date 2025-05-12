/**
# Bug in the dump-restore function

This bug occur in axisymetric case. If we try to restore the simulation, we get
a floating point exception. This doesn't occur if we don't have the library 
two-phase.h. It also doesn't occur if we don't include the library "axi.h" */ 

#include "axi.h" 
#include "navier-stokes/centered.h" 
#include "two-phase.h"

int main () {

  /**
  we are just creating a grid for the simulation, no more.*/

  init_grid(128);
  run();
}

/**
We initialise the simulation, without any specific geometry. If we had dump
a file, we are restoring the dump file. In the other case, we are just 
outputing a string in the log file.*/

event init(t=0){
  if (!restore (file = "dump"))
    fprintf(stderr, "No restoration\n");

  /**
  If we uncomment the two following line, we solve the floating point
  exception link to the restoration of the axisymetric case. */

  //  else
  //    boundary((scalar *){fm});
}

/**
We dump the simulation*/

event snapshot(i=1; i+=10; i<=100) {
  dump (file = "dump");
}

/**
#Note

This bug is not appening if we want to dump and restore the simulation in
main(), even with all the include library */
