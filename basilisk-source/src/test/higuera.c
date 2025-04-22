/**
# Viscous hydraulic Jump

We want to reproduce the hydraulic jump of Higuera (1994). */

#include "grid/cartesian1D.h"
#include "saint-venant.h"

int main() {
  X0 = 0.14;
  L0 = 1. - X0;
  G  = 1.;
  N  = 128;
  nl = 15;
  nu = 1.;
  run();
}

/**
We impose boundary condition for $h$ and $\eta$. */

h[left] = dirichlet (.2);
eta[left] = dirichlet (.2);

h[right] = dirichlet (0.02);
eta[right] = dirichlet (0.02);

/**
## Initialization

We define a scalar field *hc* to check for convergence on $h$. */

scalar hc[];

event init (i = 0) {

  /**
  We set a constant velocity at the inlet and a free outlet. */

  for (vector u in ul) {
    u.n[left] = 5.5;
    u.n[right] = neumann(0.);
  }
  
  /**
  We initialize *h* and *hc*. */
  
  foreach()
    hc[] = h[] = 0.2;
}

/**
We check for convergence. */

event logfile (t += 0.1; i <= 10000) {
  double dh = change (h, hc);
  printf ("%g %g\n", t, dh);
  if (i > 0 && dh < 1e-5)
    return 1;
}

/**
## Output

We print the elevation and the stress. */

event output (t = end) {
  vector u0 = ul[0];
  foreach()
    fprintf (stderr, "%g %g %g\n", x, eta[], 2.*u0.x[]/(h[]/nl));
}

/**
## Results

~~~gnuplot Comparison with Figure 2 of Higuera (1994).
X0 = 169
X1 = 604
Y0 = 222.24
Y1 = 528
unset tics
plot [0:][0:605] 'higuera.png' binary filetype=png with rgbimage not, \
  'log' u (X0+$1*(X1-X0)):($2/2*(Y1-Y0)+Y0) t 'h' w l,		      \
  '' u (X0+$1*(X1-X0)):($5/15*(Y1-Y0)+Y0) t 'tau' w l
~~~

## Bibliography

* Higuera, F. 1994. [The hydraulic jump in a viscous laminar
flow](http://dx.doi.org/10.1017/S0022112094002041). J. Fluid
Mech. 274, 69–92.
*/
