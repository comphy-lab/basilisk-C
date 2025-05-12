/**
# Kelvin Helmholtz 

Similar to [here](kh2.c), we setup a stably-stratified shear layer,
specified with $Ri = \frac{b_0d_0}{U_0}$. On this page, we study the
evolution with a finite disturbance, where the initial shear layer is
located at height $y_{s} = y(x) = 0 + A\mathrm{sin}(\frac{2m\pi
x}{L_0}$ .

Note that the wavelength of the distance is given by $\lambda =
\frac{L_0}{m}$. In this setup $L_0 = 100d_0$, with $d_0$ the inversion
length scale.

![KH billows emerge altough $Ri = 0.3$](kh3/b.mp4)
*/
#include "navier-stokes/centered.h"
#include "tracer.h"

scalar b[], * tracers = {b};
face vector av[];

/**
With normalized velocity and buoyancy inversions, the $d_0$ parameter
directly controls the Richardson number. 
 */

double U0 = 1, b0 = 1, d0;
int maxlevel = 9;

/**
We can specify the disturbance manually:
 */

double Ri = 0.3;
double A = 2.;
int m = 3;

/**
Or pass optional arguments to the executable in the shell:

~~~bash
$ ./a.out 3 4 0.5
~~~

for $A = 3d_0$, $\lambda = 100*\frac{d_0}{4}$ and $Ri = 0.5$.  
 */
int main(int argc, char ** argv) {
  // Read input
  if (argc >= 2)
    A = atof (argv[1]);
  if (argc >= 3) 
    m = atoi (argv[2]);
  if (argc >= 4)
    Ri = atof (argv[3]);
  assert (m > 0); // This is not allowed
  
  d0 = Ri;
  A *= d0; 
  L0 = 100*d0;
  Y0 = -L0/2;
  periodic (left);
  a = av;
  
  run();
}

/**
   The mesh in pre-refined near the shear layer before it is initialized.

   We define the distorted vertical coordinate:
   $$y' = y + A\mathrm{sin} (\frac{2m\pi x}{L_0})$$
 */

#define YACC (y + A*sin(2*pi*m*x/L0))
event init(t = 0) {
  refine (fabs(YACC) < 5*d0 && level < maxlevel - 1);
  refine (fabs(YACC) < 2*d0 && level < maxlevel);
  foreach() {
    u.x[] = U0*erf(YACC/d0);
    b[] = b0*erf(YACC/d0);
    u.y[] = U0/50.*noise()*exp(-sq(YACC/d0));
  }
  b.gradient = minmod2;
  
  /**
     For the efficient retrieval of the hydrostatic balance on a tree grid,
     the properties of the iterative solver need to tuned.
   */
  
  p.prolongation = refine_linear;
  TOLERANCE = 1e-5;
}
/**
The acceleration term due to the buoyancy field is added in the
`acceleration` event.
 */

event acceleration (i++) {
  foreach_face(y) {
    av.y[] = (b[] + b[0, -1])/2.;
  }
}
/**
Movies are generated of the buoyancy field, the adaptive mesh-refinement and the vorticity field.
 */
double tend = 35;
event movie (t += .25; t <= tend) {
  scalar lev[];
  foreach()
    lev[] = level;
  output_ppm (b, file = "b.mp4", n = 300,
	      min = -1.1, max = 1.1, linear = true);
  output_ppm (lev, file = "level.mp4", n = 300,
	      min = maxlevel - 3, max = maxlevel);
  vorticity (u, lev);
   output_ppm (lev, file = "omg.mp4", n = 300,
	       min = -U0/d0, max = U0/d0, map = cool_warm);
}


/**
The progress of the simulation is printed in the event below:
*/
event logger (i += 5) {
  printf ("%d %g %d %d %ld\n", i, t,
	  mgp.i, mgpf.i, grid->tn);
}
/**
The mesh is adapted to track the non-trivial spatial variation in the
buoyancy and the velocity-component fields.

![Evolution of the refinement level](kh2/level.mp4)
 */

event adapt (i++) {
  if (i == 10)
    TOLERANCE = 1e-3;
  double frac = 20;
  double ue = U0/frac;
  double be = b0/frac;
  adapt_wavelet ({b, u}, (double[]){be, ue, ue, ue},
		 maxlevel, maxlevel - 3);
}
