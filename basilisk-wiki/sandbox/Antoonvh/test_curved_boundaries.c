/**
# A test for curved-boundary implementations 

On this page a comparison between curved-boundary-implementations is
presented. Considering mask, Stephane's trick and the embedded
boundary method. The test case considers a two-dimensional
vortex-cylinder collision, and we consider the evolution of the
enstrophy as a critical statistic for the representation of the flow.
*/
#include "embed.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/double-projection.h"
#define RAD (pow(pow((x - xo), 2)+(pow((y - yo), 2)), 0.5))
#define ST (-(x - xo)/RAD)

#define CYLINDER (1. - sqrt(sq(x - xo) + sq(y - 5.)))
double yo = 10., xo = 12.1;
int j;

face vector muc[];
int main() {
  L0 = 25;   
  mu = muc;
  /**
     For each run, the vortex structure is placed 5 units away from
     the cylinder center. The `mask`ed run has label `j = 0`.
  */
  j = 0;
  run();
  /**
     Next, we try "Stephane's trick". 
  */
  j++;
  run();
  /**
Finally, there is the `embed`ded method.
   */
  j++;
  run();
}
/**
##Implementation

for the mask method, we need to declare a boundary internal domain
(`bid`);
*/
bid cylinder;
u.t[cylinder] = dirichlet (0.);

/**
We initialize a vortex dipole with centered location `{xo, yo}` 
*/

event init (t = 0) {
  scalar psi[];
  double k = 3.83170597;
  refine (RAD < 2.0 && level <= 9);
  refine (RAD < 1.0 && level <= 10);
  refine (fabs(CYLINDER) < 0.2 && level < 9);
  refine (fabs(CYLINDER) < 0.1 && level <= 10);
  foreach() 
    psi[] = ((RAD > 1)*ST/RAD +
	     (RAD < 1)*(-2*j1(k*RAD)*ST/(k*j0(k)) +
			RAD*ST));
  boundary ({psi});
  foreach() {
    u.x[] = -((psi[0, 1] - psi[0, -1])/(2*Delta));
    u.y[] = (psi[1, 0] - psi[-1, 0])/(2*Delta);
  }
  /**
     For the mask method, we use the `mask()' function to mask
     cells. Note that we need to unrefine the grid (i.e. implement a
     lego boundary) to get it to run.
  */
  if (j == 0) {
    unrefine (y < 7.5 && level > 7);
    mask (CYLINDER > 0. ? cylinder : none);
  }
/**
For the embedded boundary, we compute its location and implement it
like so:
 */
    if (j == 2) {
    scalar phi[];
    foreach_vertex()
      phi[] = -CYLINDER;
    fractions (phi, cs, fs);
    u.n[embed] = dirichlet (0.);
    u.t[embed] = dirichlet (0.);
  }
  boundary (all);
}
/**
We use a constant viscosity in the flow domain. This event is
compatible with all methods.
 */
event properties (i++) {
  foreach_face()
    muc.x[] = fm.x[]/500.;
  boundary ((scalar*){muc});
}
/**
   Stephane's trick is implemented via an additional event. 
 */
event Stephanes_trick (i++) {
  if (j == 1) {
    scalar f[];
    fraction (f, CYLINDER);
    foreach(){
      foreach_dimension()
	u.x[] -= u.x[]*f[];
    }
  }
}
/**
Because of the spatio-temporal localization of our problem, grid adaptation is employed.
 */
event adapt (i++)
  adapt_wavelet ({cs, u.x, u.y}, (double[]){0.001, 0.01, 0.01}, 10);
/**
##Output

First, movies are generated. 
*/
event movie ( t += 0.1 ; t <= 10){
  scalar omega[];
  vorticity (u, omega);
  foreach() {
    if (x > xo)
      omega[] = level - 5;
  }
  output_ppm (omega, n = 512, file = "movie_cyl.mp4", min = -5, max = 5);
}
/**
The overall dynamics appear quite similar. 

![The movie cycles over the three methods](test_curved_boundaries/movie_cyl.mp4)

Second, we quantify the total vorticity via the enstrophy (`E`), 
 */

event diag (i += 5) {
  double E = 0;
  boundary ({u.x, u.y});
  scalar omega[];
  vorticity (u , omega);
  foreach(){
    double vort = omega[];
    double area = dv();
    if (cs[] < 1. && cs[] > 0){ //Embedded boundary cell
      coord b, n;
      area *= embed_geometry (point, &b, &n);
      vort = embed_vorticity (point, u, b, n);
    }
    E += area*sq(vort);
  }
  char fname[99];
  sprintf (fname, "data_cyl%d", j);
  static FILE * fp = fopen (fname, "w");
  fprintf (fp, "%d\t%g\t%g\n", i, t, E);
  fflush (fp);
}
/**
~~~gnuplot The maximum value is sensitive to the boundary-implementation details, the value at t = 10 appears more robust.
set xlabel 'time'
set ylabel 'Enstrophy'
set key top right
plot 'data_cyl0' u 2:3 w l lw 3 t 'Mask',	\
     'data_cyl1' u 2:3 w l lw 3 t 'Stephan`s trick',   \
     'data_cyl2' u 2:3 w l lw 3 t 'Embedded boundary'   
~~~

   Finally we compare how long it takes for each run to complete.
 */

event stop (t = 10) {
  static FILE * fp = fopen ("perf_cyl", "w");
  timing s = timer_timing (perf.gt, iter, perf.tnc, NULL);
  fprintf (fp, "%d\t%g\t%d\t%g\n", j, s.real, i, s.speed);
  fflush (fp);
  return 1;
}

/**
Here are the performance results:

~~~gnuplot
reset
set xr [-0.5:2.5]
set xlabel 'j-label'
set ylabel 'run time [s]'
set key off
plot 'perf_cyl' u 1:2 pt 5 ps 5 
~~~

## See also

* [A comparison of four straight boundary implementations](test_straight_boundaries.c)
*/
