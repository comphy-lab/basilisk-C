/**
Running successive (identical) runs give slightly different results. It appears that the chosen timesteps differ from run to run, even if `dt` is reset. 

**EDIT:** adding the line

~~~c
if (t < 1e-10) previous = 0.;
~~~

in [timestep.h](/src/timestep.h) makes the timesteps identical in each run. The diagnostics on the field still differ from run to run (actually, they appear to differ between the first and the second run; then a ``steady state'' is reached).
*/

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

int LEVEL = 6;
int nrun = 1.;

int main()
{
  f.sigma = 1.;
  size (2.);
  origin (-1., -1.);
  for (int counter = 1; counter <= 5; counter++) {
    init_grid (1 << LEVEL);
    run();
    nrun++;
  }
}

event init (t = 0.) {
  fraction (f, -(sq(x) + sq(y/0.25) - sq(0.5)));
}

event logfile (i++) {
 if (i==1 || i ==2 || i==26) fprintf (stderr, "%d %g %g %g\n", i, t, dt, statsf(u.x).max);
}

event end (t = 0.05) {
  scalar pos[];
  foreach() pos[] = 0.;
  position (f, pos, {0,1}, add = false);
  boundary ({pos});
  double maxrad = statsf(pos).max;
  fprintf (stderr, "Run #%i ended | max =  %g\n", nrun, maxrad);
}

/**
~~~bash
$ cat test_multiple_runs/log
1 0.00028334 0.00028334 0.0384182
2 0.000823738 0.000540398 0.125714
26 0.05 0.0025419 1.84152
Run #1 ended | max =  0.126782
1 0.00028334 0.00028334 0.0364972
2 0.000823738 0.000540398 0.125097
26 0.05 0.0025419 1.84701
Run #2 ended | max =  0.126791
1 0.00028334 0.00028334 0.0364972
2 0.000823738 0.000540398 0.125097
26 0.05 0.0025419 1.84701
Run #3 ended | max =  0.126791
1 0.00028334 0.00028334 0.0364972
2 0.000823738 0.000540398 0.125097
26 0.05 0.0025419 1.84701
Run #4 ended | max =  0.126791
1 0.00028334 0.00028334 0.0364972
2 0.000823738 0.000540398 0.125097
26 0.05 0.0025419 1.84701
Run #5 ended | max =  0.126791
~~~
*/