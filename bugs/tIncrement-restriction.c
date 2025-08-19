#include "grid/multigrid.h"
#include "navier-stokes/centered.h"

// Must be compile-time constant for event system
#ifndef TFAKE
#define TFAKE 1.2e-7
#endif

double tfake = TFAKE;

int main(){  
  run();
}

event whatever (t = 0; t += tfake; t < tfake*1e2){
  fprintf(ferr, "TEST EVENT %d t=%g tfake=%g\n", i, t, tfake);
}