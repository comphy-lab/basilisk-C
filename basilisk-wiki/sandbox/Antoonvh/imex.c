/**
# Text the IMEX integrator

We solve:

$$s' = \frac{-s}{2},$$

but split it:

$$s' = \frac{s}{2} - s,$$

so that we can apply the explicit method to the first term on the RHS
and an implicit method for the second one.

~~~gnuplot 3rd order assymptotic convergence
set logscale xy
set xr [1 : 5000]
plot 'out', 3*x**(-3)
~~~
 */
#define IMEX_333 1

#include "IMEXrk.h"
scalar s[], *sl = {s};

double tend = 10;

void explicit_tendency (scalar s, scalar ds) {
  foreach()
    ds[] = s[]/2.;
}
/**
## Implicit solver:

$$s_{new} = s_{old} - dt A_{i,i} s_{new} $$

$$s_{new} = \frac{s_{old}}{1 + dt A_{i,i}}$$

This function is expected to update ($s_{old} \leftarrow s_{new}$)
aswell.
*/
void implicit_tendency (scalar s, scalar ds, double dtAii) {
  scalar si[];
  foreach() {
    si[] = s[];
    s[] = si[] / (1. + dtAii); 
  }
  foreach() 
    ds[] = (s[] - si[])/dtAii; // Generic expression 
}

int main() {
  N = 1;
  for (DT = 5e-3; DT <= 10; DT *= 2)
    run();
}

event init (t = 0) {
  foreach()
    s[] = 1;
}

event timestep (i++) {
  dt = dtnext (DT);
}

event error_eval (t = tend) {
  foreach() 
    printf ("%d %g\n", i, fabs(s[] - exp(-0.5*(t))));
  return 1;
}
