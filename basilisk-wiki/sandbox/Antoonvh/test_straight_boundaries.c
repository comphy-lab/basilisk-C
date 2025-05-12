/**
# A straight boundary in Baslisk

On this page a comparison between straight-boundary-geometry implementations is
presented. Considering the default (square domain), mask, Stephane's
trick and the embedded boundary methods. The test case considers a
two-dimensional vortex-wall collision, and we take the evolution
of the enstrophy as a critical statistic for the representation of the flow.
*/
#include "embed.h"
#include "navier-stokes/centered.h"

#define RAD (pow(pow((x - xo), 2)+(pow((y - yo), 2)), 0.5))
#define ST (-(x - xo)/RAD)
double yo, xo = 12.1;
int j;
/**
For the default and mask method, the boundary condition is implemented
at the `bottom` boundary.
 */
u.t[bottom] = dirichlet (0.);

face vector muc[];
int main() {
  L0 = 25;   
  mu = muc;
  /**
     For each run, the vortex structure is placed 5 length units away from
     the bottom. The default run has label `j = 0`. Furthermore, momentum
     conservative refinement is used.
   */
  foreach_dimension()
    u.x.refine = refine_linear;
  j = 0;
  yo = 5;
  run();
  /**
The masked boundary is implemented at a quarter of the domain, and
hence the vortex strucutre is placed higher to mimic the default run.
   */
  j++;
  yo += L0/4;
  run();
  /**
     Next, we try "Stephane's trick, where we use a slight offset. 
   */
  j ++;
  yo -= 0.001;
  run();
  /**
Finally, there is the `embed`ded method.
   */
  j++;
  run();
}
/**
## Implementation

We initialize a vortex dipole with centered location `{xo, yo}`. 
*/

event init (t = 0) {
 scalar psi[];
 double k = 3.83170597;
  refine (RAD < 2.0 && level <= 9);
  refine (RAD < 1.0 && level <= 10);
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
For the mask method, we use the `mask()` function to raise the `bottom`
boundary.
 */
  if (j == 1)
    mask (y < L0/4 ? bottom : none);
/**
For the `embed`ded boundary, we compute its location and implement it
like so:
 */
  if (j == 3) {
    scalar phi[];
    foreach_vertex()
      phi[] = y - L0/4. + 0.001;
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
  if (j == 2) {
    scalar f[];
    fraction (f, L0/4. - 0.001 - y);
    foreach(){
      foreach_dimension()
	u.x[] -= u.x[]*f[];
    }
  }
}

/**
Due to the spatio-temporal localization of our problem, grid adaptation is employed.
 */
event adapt (i++)
  adapt_wavelet ({u.x, u.y}, (double[]){0.01, 0.01}, 10);
/**
## Output

First, movies are generated, display the vorticity and the grid structure. 
*/
event movie ( t += 0.1 ; t <= 10){
  scalar omega[];
  vorticity (u, omega);
  foreach() {
    if (x > xo)
      omega[] = level - 5.;
  }
  output_ppm (omega, n = 512, file = "movie.mp4", min = -5.5, max = 5.5);
}
/**
The dynamics appear very similar.

![The movie cycles over the four methods](test_straight_boundaries/movie.mp4)

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
  sprintf (fname, "data%d", j);
  static FILE * fp = fopen (fname, "w");
  fprintf (fp, "%d\t%g\t%g\n", i, t, E);
  fflush (fp);
}
/**
~~~gnuplot The maximum value is sensitive to the method, the value at t = 10 is not.
set xlabel 'time'
set ylabel 'Enstrophy'
set key top left
plot 'data0' u 2:3 w l lw 4 t 'Square domain',	\
     'data1' u 2:3 w l lw 3 t 'Mask',         	\
     'data2' u 2:3 w l lw 2 t 'Stephane`s trick',    \
     'data3' u 2:3 w l lw 2 t 'Embedded boundary'
~~~

Finally we compare how long it takes for each run to complete.
 */

event stop (t = 10) {
  static FILE * fp = fopen("perf", "w");
  timing s = timer_timing (perf.gt, iter, perf.tnc, NULL);
  fprintf (fp, "%d\t%g\t%d\t%g\n", j, s.real, i, s.speed);
  fflush (fp);
  return 1;
}

/**
Here are the performance results:

~~~gnuplot Increasing complexity appears increasingly costsly?
reset
set xr [-0.5:3.5]
set xlabel 'j-label'
set ylabel 'run time [s]'
set key off
plot 'perf' u 1:2 pt 5 ps 3 
~~~
*/
