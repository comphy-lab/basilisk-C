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

/**
## Root cause and proposed fix

In `update_perf()` in [/src/utils.h](/src/utils.h) (lines 32–39):

~~~c
void update_perf() {
  perf.nc += grid->n;
  perf.tnc += grid->tn;
  perf.t = timer_elapsed (perf.gt);
  perf.speed = perf.tnc/perf.t;
  perf.ispeed = grid->tn/(perf.t - perf.tp);
  perf.tp = perf.t;
}
~~~

Because `timer_elapsed()` uses `gettimeofday()` (1 µs resolution) and
Basilisk enables `FE_DIVBYZERO` via `catch_fpe()`, two consecutive
iterations within the same microsecond give `perf.t == perf.tp`,
causing `grid->tn/(perf.t - perf.tp)` to raise `SIGFPE` (and similarly
`perf.t == 0.` on a very fast first iteration).

Guarding both divisions — matching the `-1.` sentinel already used by
`timer_timing()` in [/src/utils.h](/src/utils.h) (line 101) — fixes
the crash:

~~~diff
--- a/src/utils.h
+++ b/src/utils.h
@@ -33,8 +33,8 @@
   perf.nc += grid->n;
   perf.tnc += grid->tn;
   perf.t = timer_elapsed (perf.gt);
-  perf.speed = perf.tnc/perf.t;
-  perf.ispeed = grid->tn/(perf.t - perf.tp);
+  perf.speed = perf.t > 0. ? perf.tnc/perf.t : -1.;
+  perf.ispeed = perf.t > perf.tp ? grid->tn/(perf.t - perf.tp) : -1.;
   perf.tp = perf.t;
 }
~~~
*/
