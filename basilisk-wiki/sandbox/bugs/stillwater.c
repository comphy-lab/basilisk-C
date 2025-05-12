/**
This is just a simple example to demonstrate why I get the error I
reported in the forum "Problem with interaction..." (See Case#2 &
Case#3). Consider the following stationary case: When I refine the
region inside a cylinder at every single timestep using Saint-Venant.h
you can see that $h$ is not conserved anymore which is caused by
refining the cylinderical region, however $h$ is conserved exactly
when I use a uniform mesh. This does not create a big problem when I
am using Saint-Venant equations and of course when I just refine at
the initial timestep and I have a solitary wave entering the
domain. But I think this is the reason for the error I reported before
in the forum using Green-Naghdi equations since small errors in the
value of $h$ can cause large error in calculating dispersive terms
(See. Case 2 and Case 3: Case two I use adaptive mesh with minlevel=6
and maxlevel=10 and I refine the region around the cylinder at $t=0$
and the code stops at the very early timesteps due to the error in h
and u caused by refining the cylindrical regin at time step 0 while
the code works perfectly when I use $minlevel=maxlevel=10$). 

*SP: I can't see any problem. Either the problem is fixed or you need
to make the error clearer i.e. make the code crash or make a graph
showing an obvious problem.*
*/

#include "saint-venant.h"

int MINLEVEL = 6;
int MAXLEVEL = 9;

int main() {
  size (50.);
  G = 9.81;
  origin (-L0/2., -L0/2.);
  init_grid (1 << MINLEVEL);
  run();
}

double H0 = 1.0;

event init (i = 0) {
  refine (sq(x ) + sq(y)< sq(7) && level < MAXLEVEL);
  foreach() {
    h[] = H0; 
    u.x[] = 0.0; 
  }
  conserve_elevation();
}

event logfile (i++) {
  stats s = statsf (h);
  if (i == 0)
    fprintf (stderr, "t i h.min h.max h.sum u.x.rms u.x.max dt\n");
  fprintf (stderr, "%g %g %g %g %g\n", t, s.min, s.max, s.sum, dt);
}

event stop (t = 15.) {
  return 1;
}

event adapt (i++) {
  astats s = adapt_wavelet ({h}, (double[]){1e-4}, MAXLEVEL, MINLEVEL);
  refine (sq(x ) + sq(y)< sq(7) && level < MAXLEVEL);
  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
