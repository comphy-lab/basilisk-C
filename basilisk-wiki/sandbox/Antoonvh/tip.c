/**
# Solve $\ddot{x} = -x$ using inertial particles

when $x(t = 0) = 1$, the solution is the cosine function:

~~~gnuplot Looks OK for `DT` = 0.1 
set xlabel 't'
set ylabel 'x'
plot 'out' t 'discrete solution', cos(x) t 'cosine'
~~~

~~~gnuplot Convergence is also OK (first step is 1st order)
set xr [40:1000]
set logscale xy
set xlabel '# Time Steps'
set ylabel 'L_1'
plot 'log' t 'err', 100*x**-2 t '2nd order convergence'
~~~
 */
#include "inertial-particles.h"

coord p_acc (struct PA_inp inp) {
  return (coord){-inp.p.x};
}

int main() {
  L0 = 3;
  X0 = -L0/2;
  for (DT = 0.2; DT > 0.01; DT /= 2)
    run();
}

event init (t = 0) {
  new_inertial_particles (1);
  foreach_particle()
    p().x = 1;
}

event set_dtmax (i++)
  dt = dtnext (DT);

event trac (t += 0.2) {
  if (DT == 0.1) 
    foreach_particle()
      printf ("%g %g\n", t, p().x);
}

event stop (t = 10) {
  foreach_particle()
    fprintf (stderr, "%d %g\n", i, fabs(p().x - cos(t)));
}
