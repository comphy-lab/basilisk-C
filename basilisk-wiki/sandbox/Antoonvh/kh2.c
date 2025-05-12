/**
# Kelvin Helmholtz

Shearlayers are unstable to small perturbations. However, this
Kelvin-Helmholtz instability maybe surpressed by a density
stratification (and viscous effects, which are largely ignored here).

The inversion strength maybe characterized by a difference in buoyancy
$b_0$, and the shear strength by a velocity difference $U_0$. Another
important aspect of this layer is its thickness ($d_0$) as it controls
the gradients. A dimensionless number that compares the shear to the
buoyancy effects is the bulk Richardson number,

$$Ri = \frac{b_0 d_0}{U_0^2}.$$

Here we model such a stabily stratified shearlayer with $Ri = 0.1$ and $Ri = 0.4$

![Evolution of the buoyancy field for $Ri = 0.1$ and $Ri = 0.4$](kh2/b.mp4)

In a more quantitaive analysis we plot the maximum absolute vertical velocity to "detect" the instability. 

~~~gnuplot
set grid
set xlabel 'time'
set ylabel 'u_{y,max}'
plot 'out' w l lw 2 t 'Ri = 0.1',\
'log' w l lw 2 t 'Ri = 0.4'
~~~
*/
#include "navier-stokes/centered.h"
#include "tracer.h"

scalar b[], * tracers = {b};
face vector av[];

double U0 = 1, b0 = 1, d0;

int maxlevel = 8;

int main() {
  periodic (left);
   a = av;
  /**
A small value for the fluid's viscosity is set.
   */
  const face vector nu[] = {1./5000., 1./5000.};
  mu = nu;
  /**
We run the model with $Ri = 0.1$ and $Ri = 0.4$ by varying $d_0$ and
scaling the domain size accordingly. These values correspond to unstable
and stable conditions, respectively.
   */
  
  d0 = 0.1;
  L0 = 100*d0;
  Y0 = -L0/2;
  run();
  
  d0 = 0.4;
  L0 = 100*d0;
  Y0 = -L0/2;
  run();
}

/**
   The mesh in pre-refined near the shear layer before it is initialized.
 */
event init(t = 0) {
  refine (fabs(y) < 5*d0 && level < maxlevel - 1);
  refine (fabs(y) < 2*d0 && level < maxlevel);
  foreach() {
    u.x[] = U0*erf(y/d0);
    b[] = b0*erf(y/d0);
    u.y[] = U0/100.*noise()*exp(-sq(y/d0));
  }
  /**
For the efficient retrieval of the hydrostatic balance on a tree grid,
the properties of the iterative solver need to tuned.
   */
  p.prolongation = refine_linear;
  TOLERANCE = 1e-5;
}
/**
The acceleration term of the buoyancy field is added in the
`acceleration` event.
 */

event acceleration (i++) {
  foreach_face(y) {
    av.y[] = (b[] + b[0, -1])/2.;
  }
}
/**
Movies are generated of the buoyancy field and the adaptive mesh-refinement.
 */
double tend = 30;
event movie (t += .25; t <= tend) {
  scalar lev[];
  foreach()
    lev[] = level;
  output_ppm (b, file = "b.mp4", n = 300,
	      min = -1., max = 1., linear = true);
  output_ppm (lev, file = "level.mp4", n = 300,
	      min = maxlevel - 3, max = maxlevel);
}

/**
The maximum of the absolute vertical velocity is logged to trace the
evolution of the instability.
 */
event KH_tracing (i++) {
  FILE * fp;
  if (d0*b0/sq(U0) > 0.25)
    fp = stdout;
  else
    fp = stderr;
  fprintf (fp, "%g %g\n", t, normf(u.y).max);
}
/**
The mesh is adapted to track the features of the buoyancy and the
velocity-component fields.

![Evolution of the refinement level](kh2/level.mp4)
 */

event adapt (i++) {
  double frac = 20;
  double ue = U0/frac;
  double be = b0/frac;
  adapt_wavelet ({b, u}, (double[]){be, ue, ue, ue},
		 maxlevel, maxlevel - 3);
}
