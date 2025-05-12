/**
# The [CBL case](http://www.basilisk.fr/sandbox/Antoonvh/vhm16.c) for testing SGS TKE scheme from Antoon. 

*/
#include "grid/octree.h" // <- Uncomment for 3D
#include "navier-stokes/centered.h"
#include "../Antoonvh/profile5b.h" // A grid-adaptive profiling function
#define damp(s) (-(s) * (exp((y - 2)) - 1) * (y > 2) / 2.)
#define Emin 5E-5

scalar b[];
// scalar * tracers = {b};
scalar e120[];
scalar *tracers = {b, e120};
#include "SGS_TKE.h"
face vector av[];

double ue = 0.02; //U / ?
double be = 0.02; //b_0 / 50
int maxlevel = 6;

/**
We use normalized values for the initital stratification strength
($\mathrm{d}b/\mathrm{d}y = N^2$) and the bottom buoyancy ($b_0$). As
such that the CBL-height lengthscale $\mathcal{L} = 1$. The Prandtl
number ($Pr$) also has a value of unity and the Reynolds number ($Re$)
is defined below as `Re`$=3000$.
*/
double Re = 3000;

double nu, T = 1, Tend = 1, dT = 1;

u.t[bottom] = dirichlet(0.);
b[bottom]   = dirichlet(1.);
b[top]      = neumann(1.);
int main() {
  periodic (left);
  nu = pow(1./Re, 9./8.);
#if (dimension > 2)
  Evis[bottom] = dirichlet(0.); // Flux is explicitly calculated
  Evis[top] = dirichlet(0.);
  periodic (front);
  u.r[bottom] = dirichlet(0.);
#endif
  T  = pow(nu, -1./3.);
  dT = T;
  Tend = 30 * T;
  L0 = 3.;
  a = av;
  /**
     The most important updates are related to using non-default
     attributes for the prognostic variable fields.
  */
  foreach_dimension()
    u.x.refine = refine_linear; // Momentum conservative refinement
  p.refine = p.prolongation = refine_linear; // Third order accurate for vertical variations
  b.gradient = minmod2; // A flux limiter
  init_grid (1 << 5);
  run();
}

event init (t = 0) {
  TOLERANCE = 10E-5;
  DT = 0.01; 
  refine (y < 0.1  && level < (maxlevel - 1));
  refine (y < 0.05 && level < maxlevel);
  foreach(){
    b[] = y + 0.001 * noise();
    e120[] = (y < 1) ? 0.1 : Emin;
  }
}

event acceleration (i++) { //Gavity and damping for u.y
  foreach_face(y)
    av.y[] = (b[] + b[0,-1])/2. + damp((u.y[] + u.y[-1])/2.);
}

event tracer_diffusion (i++) { // and damping for b.
  diffusion (b, dt, Kh);
  foreach()
    b[] += damp(b[] - y)*dt;
}

event adapt (i++) {
  DT = min(DT * 1.05, 0.1); // DT_max = N/10
  TOLERANCE = min(TOLERANCE * 1.05, 10E-3);
  adapt_wavelet((scalar *){b, u}, (double[]){be, ue, ue, ue}, maxlevel);
}

/**
The keep the required computational resources within limited bounds,
the simulation is stopped at $t=$`Tend`.
 */
event stop (t = Tend){
  return 1;
}

/**
## Output

The output consists of vertical profiles,
*/
event profiler (t += dT) {
  char fname[99];
  sprintf (fname, "proft=%dT", (int)(t/T));
  profile ((scalar*){b, u}, fname);
}

event progress(i++)
{
  fprintf(stderr, "i=%d t=%g p=%d u=%d \n", i, t, mgp.i, mgu.i);
}

/**
time series data for some domain-integrated quantities and solver
characteristics,
 */
event time_series (i += 25) {
  double e = 0, diss = 0, sbflx = 0;
  foreach(reduction(+:diss) reduction(+:e) reduction(+:sbflx)) {
    foreach_dimension() {
      e += sq(u.x[])*dv();
      diss += dv()*(sq(u.x[1] - u.x[-1]) +
		    sq(u.x[0,1] - u.x[0,-1]) +
		    sq(u.x[0,0,1] - u.x[0,0,-1]))/sq(2.*Delta);
    }
    if (Delta > y)
      sbflx += (b[0,-1] - b[])*dv()/sq(Delta);
  }
  diss *= -nu;
  sbflx *= nu;
  e /= 2.;
  static FILE * fp = fopen ("timeseries", "w");
  if (i == 0)
    fprintf (fp, "t\ti\tn\twct\tspeed\te\tdiss\tsbflx\n");
  fprintf (fp, "%g\t%d\t%ld\t%g\t%g\t%g\t%g\t%g\n",
	   t/T, i, grid->tn, perf.t, perf.speed, e, diss, sbflx);
  fflush (fp);
}

/**
   Furthermore, two movies are rendered that display the evolution of
   the buoyancy structures (using the field `m` $=\mathrm{ln}\left( 
   \| \nabla b \|+1 \right)$), and the numerical mesh.
*/
event movies (t += 0.25) {
  vector db[];
  scalar m[];
  boundary ({b});
  gradients ({b}, {db});
  foreach() {
    m[] = 0;
    foreach_dimension()
      m[] += sq(db.x[]);
    if (m[] > 0)
      m[] = log(sqrt(m[])+1.);
  }
  boundary ({m});
  output_ppm (m, file = "db.mp4", n = (1 << maxlevel), min = 0, max = 2, linear = true);
  scalar lev[];
  foreach()
    lev[] = level;
  output_ppm (lev, file = "lev.mp4", n = (1 << maxlevel), min = 1, max = maxlevel);
}

/**
## Results

Here, the simulation in run in 2D and these are the two movies, 

![The evolution of the aforementioned `m` field](vhm16/db.mp4)

![The evolution of the grid structure](vhm16/lev.mp4)

Some profiles of the buoyancy field.

~~~gnuplot
set xr [0:2]
set yr [0:2]
set xlabel 'Buoyancy' font ",15"
set ylabel 'height' font ",15"
set key box top left
set size ratio -1
set grid
plot 'proft=1T' u 2:1 w l lw 3 t 't = 1T' , \
'proft=3T' u 2:1 w l lw 3 t 't = 3T' ,	    \
'proft=8T' u 2:1 w l lw 3 t 't = 8T' ,	    \
'proft=20T' u 2:1 w l lw 3 t 't = 20T'
~~~

The time evolution of the total kinetic energy:

~~~gnuplot
set xr [0 : 30]
set yr [0 : 0.022]
set xlabel 'time/T'
set ylabel 'Energy'
set key off
set size square
plot 'timeseries' u 1:6 w l lw 3 
~~~

It appears that the decay of 2D turbulence differs quite a bit from it's 3D
counterpart, see also [this
example](http://basilisk.fr/examples/turbulence.c).

# References

Van Heerwaarden, C. C., & Mellado, J. P. (2016). Growth and decay of a
convective boundary layer over a surface with a constant
temperature. Journal of the Atmospheric Sciences, 73(5), 2165-2177.

van Hooft, J. A., Popinet, S., van Heerwaarden, C. C., van der Linden,
S. J., de Roode, S. R., & van de Wiel, B. J. (2018). Towards Adaptive
Grids for Atmospheric Boundary-Layer Simulations. Boundary-Layer
Meteorology, 167(3), 421-443.
*/