/** 

# Memory leak in solve.h

We use the kuramoto.c example to test for memory leaks in solve.h. At each 
iteration, the memory used by realloc scalar increases. To see this, we run with 

~~~bash
CFLAGS="-O2 -DMTRACE=1" make memory_leak_solve.tst
~~~

and plot memory use as usual. We can see that `realloc_scalar()` increases at
each step. A possible fix is to add `delete ({_res, _da});` at the end of
the `solve()` macro.

~~~
diff -rN -u old-basilisk/src/solve.h new-basilisk/src/solve.h
--- old-basilisk/src/solve.h	2026-02-20 10:09:11.808248041 +0100
+++ new-basilisk/src/solve.h	2026-02-20 10:09:11.808421631 +0100
@@ -172,6 +172,7 @@
 	     "  res: %g nrelax: %d\n", LINENO, a.name,
 	     _s.i, _s.resa, _s.nrelax),
       fflush (ferr);
+  delete ({_res, _da});
   return _s;
 }}
~~~

*/



#include "grid/multigrid1D.h"
#include "solve.h"

int main()
{  
  init_grid (128);
  L0 = 32.*pi;
  periodic (right);
  scalar u[];
  foreach()
    u[] = cos(x/16.)*(1. + sin(x/16.));
  
  double dt = 1e-1;
  int i = 0;
  TOLERANCE = 1e-6;
  for (double t = 0; t <= 10; t += dt, i++) {
    if (i % 1 == 0) {
      foreach()
        fprintf (stdout, "%g %g %g\n", t, x, u[]);
      fputs ("\n", stdout);
    }
    scalar b[];
    foreach()
      b[] = u[] - dt*u[]*(u[1] - u[-1])/(2.*Delta);
    mgstats stats = solve (u,
			   u[] + dt*(u[-1] - 2.*u[] + u[1])/sq(Delta)
			   + dt*(u[-2] - 4.*u[-1] + 6.*u[] - 4.*u[1] + u[2])/sq(sq(Delta)),
			   b[]);
    fprintf (stderr, "%g %d\n", t, stats.i);
  }
}