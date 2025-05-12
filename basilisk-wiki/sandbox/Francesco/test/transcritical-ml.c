/**
# Transcritical flow over a bump with multiple layers

We want to reproduce the transcritical test case of [Audusse et al, 2011](/src/references.bib#audusse2011). The viscosity is set to $\nu = 0.01 m^2/s$ and for the bottom friction we use the Strickler relation:
$$
k(h,\mathbf{U}) = \frac{g}{S^2h^{1/3}}|\mathbf{U}| 
$$
with $S = 25 m^{1/3}/s$ the Strickler coefficient, h the water depth and $\mathbf{U}$ the depth-averaged velocity.*/ 

#include "grid/cartesian1D.h"
#include "saint-venant.h"

int main() {
  X0 = 0.;
  L0 = 21.;
  G = 9.81;
  N = 256;
  nu = 0.01;
  gradient = zero;

  /**
  We do the test for several number of layers.*/
  
  for (nl = 2; nl <= 15; nl++)
    run();
}

/**
We impose the outlet water level. */

h[right] = dirichlet(0.6);
eta[right] = dirichlet(0.6);

/**
## Initialitazion

We initialize the topography, the *h* and we create a field *hc* to check convergence on *h*. */

scalar hc[];

event init (i = 0) {
  foreach() {
    zb[] = max(0., 0.2*(1. - 1./sq(5.75/2.)*sq(x - 10.)));
    h[]  = 0.6 - zb[];
    hc[] = h[];
  }

  /**
  ## Boundary Conditions on velocity
  
  We impose a constant inflow of 1 m^2/s in the inlet and neumann condition on the outlet.*/
  
  for (vector u in ul) {
    u.n[left] = dirichlet(h[left] ? 1./h[left] : 0.);
    u.n[right] = neumann(0.);
  }
}

/** We check for convergence.*/

event logfile (t += 0.1; i <= 100000) {
  double dh = change (h, hc);
  if (i > 0 && dh < 1e-3)
    return 1;
}

/* uncomment this part if you want on-the-fly animation. */
event output (i++) {
  static FILE * fp = popen ("gnuplot", "w");
  fprintf (fp, 
           "set title 't=%f'\n"
           "set xl 'x'\nset yl 'h'\n"
           "plot [0:21][] '-' w l t 'eta', '-' u 1:3 w l t 'zb'\n", t); 
  foreach()
    fprintf (fp, "%g %g %g\n", x, eta[], zb[]);
  fprintf (fp, "e\n");
}

/**
## Output

We save the value of $\eta$.*/

event output (t = end) {
  FILE * fp1 = fopen ("end", "w");
  foreach()
    fprintf (fp1, "%g %g %g\n", x, eta[], zb[]);
  fprintf (fp1, "\n\n");
}

/**
## Results

~~~gnuplot Free surface and topography.
set xr [0:21]
set yr [0:1]
set xlabel 'x'
set ylabel 'eta'
plot 'end' i 0 w l t '2 layers', '' i 3 w l t '5 layers', '' i 13 w l t '15 layers', '' i 0 u 1:3 w l t 'topography' 
~~~
*/

