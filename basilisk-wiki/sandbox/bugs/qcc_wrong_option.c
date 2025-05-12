/**

When running this code with a non-existing option of qcc (tested using a Makefile), for example using

qcc ... -disable-dimsensions ...

where "disable-dimsensions" is wrong due to a typo, qcc hangs and doesn't return an error.
It has also to be noted that some files ( *.expand, *.vregs, *.into_clglayout, *.outof_clglayout, *.jump, *.reginfo) are created and indefinitely increase in size.

*/


#include "compressible/two-phase.h"
#include "compressible/Mie-Gruneisen.h"

int main(){
  init_grid(1<<1);
  run();
}

event end(i=1);