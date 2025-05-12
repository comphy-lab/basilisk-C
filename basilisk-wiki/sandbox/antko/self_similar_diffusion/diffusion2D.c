/**
# Self-similar diffusion in 2D

A 2D blob of concentration $c(r,t)$ initially occupying the unit disk freely diffuses according to:
$$
\frac{\partial c}{\partial t} = \Delta c.
$$
The concentration is initially constant and evolves under the constraint:
$$
2 \pi \int c(r,t) \, r \mathrm{d}r = 1
$$*/

#if 0
  #include "axi.h"
#endif
#include "run.h"
#include "diffusion.h"

scalar concentration[];
face vector diffusion_coef[];
const double tend = 10.;

mgstats mgconcentration;

/** At the box boundaries, the concentration is set to 0. */

concentration[left]   = dirichlet(0.);
concentration[right]  = dirichlet(0.);
concentration[top]    = dirichlet(0.);

/** The computational box is ten times larger than the initial patch, and the number of points is set to ensure that there are at least ten points in the initial patch*/

int main() {
  size(20.);
  origin(-10.,0.);
  N = 256;
  init_grid (N);
  run();
}

/** The initial value for the concentration is chosen to satisfy the integral constraint.*/

#define circle(x,y) (sq(1.) - sq(x) -sq(y))

event init (i = 0) {
  #if AXI
  foreach()
    concentration[] = (circle(x,y) > 0.) ? 1. / (4./3.*M_PI) : 0.;
  #else
  foreach()
    concentration[] = (circle(x,y) > 0.) ? 1. / M_PI : 0.;
  #endif
  boundary({concentration});
  
/** The diffusion coefficient is set to 1.*/
  
  foreach_face()
    diffusion_coef.x[] = fm.x[];
  boundary ((scalar *){diffusion_coef});
}

/** The maximal time step is the diffusive timestep. */

event integration (i++) {
  dt = dtnext (sq(L0/N)/2.);
  mgconcentration = diffusion (concentration, dt, D = diffusion_coef);
}

/** Regularly we output the concentration profile, in native variables and in the self-similar space.*/

event output (t+=0.5; t<tend) {
  static FILE * fp = fopen("cprof", "w");
  static FILE * fpr = fopen("cprof_resc", "w");
  scalar concentration_rescaled[];
  #if AXI
  foreach()
    concentration_rescaled[] = concentration[] * pow(t, 3./2.);
  #else
  foreach()
    concentration_rescaled[] = concentration[] * t;
  #endif
  boundary({concentration_rescaled});
  for (double y = 0.; y <= 10.; y += 0.01) {
    fprintf (fp, "%g %g\n", y,
      interpolate (concentration, 0., y));
    if (t > 0) fprintf (fpr, "%g %g\n", y / sqrt(t),
      interpolate (concentration_rescaled, 0., y));
  }
  fprintf (fp, "\n");
  fprintf (fpr, "\n");
  if (t == tend) {
    fclose (fp);
    fclose (fpr);
    static FILE * fpc = fopen ("concentration.ppm", "w");
    output_ppm (concentration, fpc);
    fclose (fpc);
  }
}

event end (t=tend) {

}

/**
## Results

The concentration profiles are reported in the physical space:

~~~gnuplot Raw concentration profiles.
plot 'cprof' w l
~~~

and in the self-similar space. The self-similar solution for the diffusion equation,
$$
c(r,t) = \frac{1}{4 \pi t} e^{-r^2/4t}
$$
is also represented for comparison:

~~~gnuplot Concentration profiles in the self-similar space.
plot 'cprof_resc' w l, 1./(4*pi)*exp(-x*x/4.)
~~~
*/