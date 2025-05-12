/**
# Transcritical flow over a bump

The water, flowing over a bump, goes from subcritical to
supercritical and then return again subcritical. The analytical 
solution was obtained thanks to the SWASHES library. */

#include "grid/cartesian1D.h"
#include "saint-venant.h"

int main() {
  L0 = 25;
  G = 9.81;
  DT=HUGE;
  for (N = 32; N <= 256; N *= 2)
    run();
}

/**
We impose the outlet water level and a constaint inflow. */
u.n[left] = dirichlet(h[left] ? 0.18/h[left] : 0.);
u.n[right] = neumann(0.);
h[right] = 0.33;
eta[right] = 0.33;

/**
## Initialisation

We initialise the topography, the water depth *h* and we
create a field *hc* to check convergence on *h*. */

scalar hc[];

event init (i = 0) {
  foreach() {
    zb[] = max(0., 0.2 - 0.05*(x - 10.)*(x - 10.));
    hc[] = h[] = 0.33 - zb[];
  }
}

/**
We check for convergence. */

event logfile (t += 0.1; i <= 100000) {
  double dh = change (h, hc);
  if (i > 0 && dh < 1e-5)
    return 1;
}

/**
## Outputs*/

event output (t = end) {
  char name[80];
  sprintf (name, "end-%d", N);
  FILE * fp = fopen (name, "w");
  foreach()
    fprintf (fp, "%g %g %g\n", x, eta[], zb[]);
  fclose(fp);
}

/**
## Results

~~~gnuplot Numerical and analytical solutions for transcritical flow over a bump. 
set xr [5:15]
set xlabel 'x'
set ylabel 'z'
plot 'swashes' u 1:4 t 'topography' w l, \
     'swashes' u 1:6 t 'free surface (analytical)' w l, \
     'end-256' w p t 'free surface (numerical)'
~~~

## See also

* [Same test with Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/shock.html)
*/
