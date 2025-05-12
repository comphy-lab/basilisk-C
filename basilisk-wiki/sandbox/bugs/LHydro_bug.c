/**
# Minimal model for bug in multi-layer model

*This is not a bug, or rather it is a bug in your code. You are overwriting the pointer 'H' to a field you allocated (in `out0`) with an already allocated field (`hl[0]`). This causes the problem which you correctly describe below. The moral of the story is that you have to be careful when manipulating field pointers / lists.*

compile with
qcc -autolink -g -lm -o LHydro_bug LHydro_bug.c
Model crashes with segmentation fault. What I think happens is that at the end of event out0 the model deletes hl[0] as well as the scalar created in the event.

Does not happen if layered/remap.h is not included
Does not happen if the 2 lines in out0 are replaced by the commented out line
*/

#include "layered/hydro.h"
#include "layered/remap.h"

int main()
{
  nl=2;
  init_grid (1 << 6);
  run();
}

event out0 (t = 0.)
{
  scalar H[];
  H = hl[0];
  //scalar H = hl[0];
}
