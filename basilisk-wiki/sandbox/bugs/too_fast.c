/**
# Basilisk cant iterate faster than the timer can tick

The speed performance computation in `common.h` can attempt a devide-by-zero, if a run loop iteration isnt very demanding. This can happen in test cases. e.g. 

[https://basilisk.fr/sandbox/Antoonvh/dagan_fig2a.c]()

or in a "do nothing" program:
*/
#include "run.h"

int main() {
  DT = 1;
  run();
}

event set_dt (i++) {
  dt = dtnext(DT);
}

event stop (t = 100);
