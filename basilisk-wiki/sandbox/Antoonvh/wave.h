/**
# A Wave-equation solver

An solver for, 

$$\frac{\partial^2 s}{\partial t^2} = \nabla^2s + a,$$

Using Newmark's method (in implicit trapezoidal mode),

 */
#include "run.h"
#include "poisson.h"

// The solution and its time derivative
scalar s[], ds[];
scalar a[];


//Implicit trapezoidal scheme coefficients for Newmark's method
double beta = 0.25, gamma1 = 0.5;


// Update ds/dt using the Laplacian of s.
void advance_ds (scalar s, scalar ds, double dt) {
  foreach() {
    double lap = 0;
    foreach_dimension()
      lap += (-2*s[] + s[1] + s[-1])/sq(Delta);
    ds[] += dt*(lap + a[]);
  }
}

event dt_setter (i++) {
  dt = dtnext(DT);
}

event acceleration (i++);

event integrate (i++) {
  scalar b[]; //rhs
  const scalar lam[] = -1;
  const face vector alpha[] = {beta*sq(dt), beta*sq(dt)};
  
  foreach() {
    double lap = 0;
    foreach_dimension()
      lap += (-2*s[] + s[1] + s[-1])/sq(Delta);
    b[] = -s[] - dt*ds[] - sq(dt)/2.*(1 - 2*beta)*lap - sq(dt)*a[]/2.;
  }
  advance_ds (s, ds, dt*(1 - gamma1));
   poisson (s, b, alpha, lam);
  advance_ds (s, ds, dt*gamma1);
}

/**
## Usage

"Usage"

* [A 1D travelling wave](wave.c)
* [A Gaussian bump](bump_wave.c)
* [A source of waves](source.c)
 */
