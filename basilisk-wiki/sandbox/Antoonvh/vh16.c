/**
# The Growth and Decay of Atmospherc Convective Turbulence.

Here we follow Van Heerwaarden and Mellado (2016) and model the growth
and subsequent decay of convective turbulence in the atmosphere. This
is very similar as was done in Van Hooft et al. (2018). However, there
have been some advancements in the set-up methodology.
*/
//#include "grid/octree.h" // <- Uncomment for 3D
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "profile5b.h"
#include "profiling.h"

#define damp(s) (-(s)*(exp((y-2))-1)*(y>2)/2.)

scalar b[];
scalar * tracers = {b};
face vector av[];

double ue = 0.02; //U / ?
double be = 0.02; //b_0 / 50
int maxlevel = 9;

double Re = 3000;

double nu, T = 1, Tend = 1, dT = 1;

u.t[bottom] = dirichlet(0.);
b[bottom] = dirichlet(1.);
b[top] = neumann(1.);
int main(){
  periodic(left);
  nu = pow(1./Re, 9./8.);
#if (dimension == 2)
  const face vector muc[] = {nu, nu};
#elif (dimension > 0)
  periodic(front);
  const face vector muc[] = {nu, nu, nu}
#endif
  T = pow(nu, -1./3.);
  dT = T;
  Tend = 30 * T;
  L0 = 3.;
  mu = muc;
  a = av;
  /**
     The most important updates are related to using non-default attributes for the prognostic variable field. 
  */
  foreach_dimension()
    u.x.refine = refine_linear;
  p.refine = p.prolongation = refine_linear;
  b.gradient = minmod2;
  init_grid (1<<5);
  run();
}

event init (t = 0){
  TOLERANCE = 10E-5;
  DT = 0.01; 
  refine(y< 0.1 && level < (maxlevel - 1));
  refine(y< 0.05 && level < maxlevel);
  foreach()
    b[] = y + 0.001*noise();
}

event acceleration(i++){ //Gavity and damping for u.y
  foreach_face(y)
    av.y[] = (b[]+b[0,-1])/2. + damp((u.y[]+u.y[-1])/2.);
}

event tracer_diffusion(i++){ // and damping for b.
  diffusion(b, dt, mu);
  foreach()
    b[] += damp(b[]-y)*dt;
}

event adapt(i++){
  DT = min(DT * 1.05, 0.1); // DT_max = N/10
  TOLERANCE = min(TOLERANCE * 1.05, 10E-3);
  adapt_wavelet((scalar *){b,u},(double[]){be, ue, ue, ue}, maxlevel);
}

/**
The keep the required computational resources within limited bounds,
the simulation is stopped at $t=$`Tend`.
 */
event stop(t = Tend){
  return 1;
}

/**
## Output

The output consists of vertical profiles,
*/
event profiler(t += dT){
  char fname[99];
  sprintf (fname, "proft=%dT", (int)(t/T));
  profile ((scalar*){b, u}, fname);
}

/**
time series data for some domain-integrated quantities and solver characteristics, 
 */
event time_series (i += 25){
  double e = 0, diss = 0, sbflx = 0;
  foreach(reduction(+:diss) reduction(+:e) reduction(+:sbflx)){
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
  static FILE * fp = fopen("timeseries", "w");
  if (i==0)
    fprintf(fp, "t\ti\tn\twct\tspeed\te\tdiss\tsbflx\n");
  fprintf(fp, "%g\t%d\t%ld\t%g\t%g\t%g\t%g\t%g\t%g\n",
	  t/T, i, grid->tn, perf.t, perf.speed, e, diss, sbflx);
  fflush(fp);
}
/**
and simulation dumps for run restoration and/or post processing.
*/

event dumping (t += T){
  char fname[99];
  sprintf(fname, "dump_t=%dT",(int)((t/T)+0.5));
  dump(fname);
}

/**
   Furthermore, two movies are rendered that display the evolution of
   the byoyancy structures (using the field `m` $=\mathrm{log}\left(
   \| \nabla b \|+1 \left)$), and the numerical mesh.
*/
#if (dimension == 2)
event movies(t += 0.25){
  char fname[99];
  vector db[];
  scalar m[];
  boundary({b});
  gradients ({b}, {db});
  foreach(){
    m[] = 0;
    foreach_dimension()
      m[] += sq(db.x[]);
    if (m[] > 0)
      m[] = log(sqrt(m[])+1.);
  }
  boundary({m});
  sprintf(fname , "ppm2mp4 db.mp4");
  static FILE * fpm = popen(fname, "w");
  output_ppm (m, fpm, n = (1<<maxlevel), min = 0, max = 2, linear = true);

  scalar lev[];
  foreach()
    lev[] = level;
  sprintf(fname , "ppm2mp4 lev.mp4");
  static FILE * fpl = popen(fname, "w");
  output_ppm (lev, fpl, n = (1<<maxlevel), min = 1, max = maxlevel);
}
#endif

/**
## Results

Here, the simulation in run in 2D and these are the two movies, 

![The evolution of the aforementioned `m` field](vh16/db.mp4)

![The evolution of the grid structure](vh16/lev.mp4)

Some profiles of the buoyancy field.

~~~gnuplot
set xr [0:2]
set yr [0:2]
set xlabel 'Buoyancy'
set ylabel 'height'
set key box top left
set size ratio -1
plot 'proft=1T' u 2:1 w l lw 3 t 't = 1T' , \
'proft=3T' u 2:1 w l lw 3 t 't = 3T' ,	    \
'proft=8T' u 2:1 w l lw 3 t 't = 8T' ,	    \
'proft=20T' u 2:1 w l lw 3 t 't = 20T'
~~~

Here is the time evolution of the total kinetic energy.

~~~gnuplot
reset xr
reset yr 
set xlabel 'time/T'
set ylabel 'Energy'
set key off
set size square
plot 'timeseries' u 1:6 w l lw 3 
~~~

It appears the decay of 2D turbulence differs quite a bit from its 3D
counterpart, see also [this
example](http://basilisk.fr/examples/turbulence.c).


## References

Van Heerwaarden, C. C., & Mellado, J. P. (2016). Growth and decay of a
convective boundary layer over a surface with a constant
temperature. Journal of the Atmospheric Sciences, 73(5), 2165-2177.

van Hooft, J. A., Popinet, S., van Heerwaarden, C. C., van der Linden,
S. J., de Roode, S. R., & van de Wiel, B. J. (2018). Towards Adaptive
Grids for Atmospheric Boundary-Layer Simulations. Boundary-Layer
Meteorology, 167(3), 421-443.
*/
