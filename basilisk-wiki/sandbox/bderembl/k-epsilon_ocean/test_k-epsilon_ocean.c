
#include "grid/multigrid.h"
#include "layered/hydro.h"
#define EOS_LINEAR 1
#include "eos_ocean.h"
#include "dr_TS.h" // proposition to have T and S in dr.h
#include "k-epsilon_ocean.h"



int main() {

  foreach_dimension()
    periodic (right);

  N = 1;
  nl = 50;

  init_grid (N);  
}

