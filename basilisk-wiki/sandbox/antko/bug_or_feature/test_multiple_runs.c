#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

int LEVEL = 6;
int nrun = 1.;

int main()
{
  f.sigma = 1.;
  size (2);
  origin (-1., -1.);
  for (int counter = 1; counter <= 5; counter++) {
    init_grid (1 << LEVEL);
    run();
    nrun++;
  }
}

event init (t = 0.) {
  fraction (f, -(sq(x) + sq(y/0.25) - sq(1.)));
}

event end (t = 1.) {
  scalar pos[];
  position (f, pos, {0,1});
  double maxrad = statsf(pos).max;
  fprintf (stderr, "Run #%i : max = %g\n", nrun, maxrad);
}

/**
We would expect the same value at the end of each run, but we get:

~~~bash
$ cat test_multiple_runs/out
Run #1 : max = 0.831405
Run #2 : max = 0.940244
Run #3 : max = 0.733233
Run #4 : max = 0.733233
Run #5 : max = 0.733233
~~~
*/
