/**
# A rising moist bubble

A moist bubble is placed in the atmosphere and rises because water
vapor is lighter than air. It reaches condensation levels soon to form
a cloud. The setup is *inspired* by Grabowski and Clark (1991) and the
[setup of B van
Stratum](https://github.com/Chiil/microhh/tree/microphysics/).
 */
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "thermo.h"
#include "bwatch.h"

int maxlevel = 8;

int main(int argc, char * argv[]) {
  if (argc > 1)
    maxlevel = atoi (argv[1]);
  const face vector K[] = {0.1, 0.1, 0.1};
  mu = K;
  L0 = 3600;
  X0 = Z0 = -L0/2;
  periodic (left);
#if (dimension == 3)
  periodic (front);
#endif
  N = 64;
  run();
}

#define RADIUS (sqrt(sq(x) + sq(y - 800.) + sq (z)))
event init (t = 0) {
  T_ref = 283.;
  P0 = 1e5;
  double dlnthetadz = 1.0e-5;
  refine (RADIUS < 350. && level < maxlevel);
  foreach() 
    thl[] = exp(log(283.) + dlnthetadz*y) + 0.05*noise();
  /**
     $q_t$ may be initialized via the relative humidity (using the
     `QSAT` macro). However it requires an iterative method to find
     the corresponding hydrostatic pressure profile. The procedure is
     started after we have taken a first guess.
  */
  scalar p_temp[];
  int j = 0, j_max = 10;
  set_pres (guess = true);
  do {
    foreach() {
      double r = RADIUS;
      if (r < 200)
	qt[] = QSAT;
      else if (r < 300) 
	qt[] = QSAT * (0.2 + 0.8 * sq(cos(pi*(r - 200.)/200.)));
      else 
	qt[] = QSAT * 0.2;
    }
    set_pres ();
    j++;
  } while (change (pres, p_temp) > .01 && j < j_max);
  if (j == j_max)
    fprintf (stderr,
	     "Static Pressure not found with %d iterations\n", j);
  system ("wget https://previews.123rf.com/images/naypong/naypong1707/naypong170700178/81549848-top-view-of-natural-green-grass-texture-aerial-view-of-park.jpg");
  system ("convert 81549848-top-view-of-natural-green-grass-texture-aerial-view-of-park.jpg -resize 400x400! grass.jpg");
}

#include "diffusion.h"
event tracer_diffusion (i++) {
  for (scalar s in tracers)
    diffusion (s, dt, mu);
}

event logger (i += 5) 
  fprintf(stdout, "%g %d\n", t, i);

event adapt (i++) 
  adapt_wavelet ((scalar*){qt, u},
		 (double[]){0.0005, 0.1, 0.1, 0.1}, maxlevel, 5);

event stop (t = 1000);
/**
## Output and Results

We output a movie. 

![](moist_bubble/cloud.mp4)
 */
event movie (t += 5) {
  scalar qc[];
  foreach() 
    qc[] = QC;
  static FILE * fp = popen ("ppm2mp4 cloud.mp4", "w");

  default_lights();
  lights[1].dir.y *= 2;
  normalize (&lights[2].dir);
  watch (fov = 2500, O = {5000, 3000, 2000},
	 poi = {0, 1000, 0});
  image ("grass.jpg", n = {0,1,0}, alpha = 1);
  volume (qc, sc = 0.02, cols = true, shading = 1, min = -1, max = 1);
  store (fp);
  plain();
}


/**
##Reference

Grabowski, Wojciech W., and Terry L. Clark. "Cloud–environment
interface instability: Rising thermal calculations in two spatial
dimensions." Journal of the Atmospheric Sciences 48.4 (1991): 527-546.
 */
