/**

qcc hangs forever instead of returning an error when the EOS is missing when using the compressible solver.
It corresponds to the commented line below.

*/


#include "compressible/two-phase.h"
//#include "compressible/Mie-Gruneisen.h"  // EOS, qcc hangs forever if missing !!

int main(){
  init_grid(1<<1);
  run();
}

event end(i=1);