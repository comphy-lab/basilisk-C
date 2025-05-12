/**
# Testing adaption */

#include "spherical.h"
#include "saint-venant.h"

/**
We then define a few useful macros and constants. */

#define MAXLEVEL 17
#define MINLEVEL 5
#define ETAE     1e-2 // error on free surface elevation (1 cm)
#define LON0  0.
#define LAT0   0.
#define DOMAIN_SIZE 1.

int main()
{

  Radius = 6371220.;
  // the domain is 1 degrees squared
  size (DOMAIN_SIZE);
  // centered on 0,0 longitude,latitude
  origin (LON0 - L0/2.,LAT0  - L0/2.);

  init_grid (1 << MINLEVEL);

  /**
  We then call the *run()* method of the Saint-Venant solver to
  perform the integration. */

  run();
}

scalar lim[];

/**
## Adaptation
*/
int adapt() {
  foreach(){
    double bound = L0*0.05;
    double xlim = fabs(x-LON0) < bound ? 1. : 0.;
    double ylim = fabs(y-LAT0) < bound ? 1. : 0.;
    lim[] = xlim * ylim ;}
  boundary ({lim});

  /**
  We can now use wavelet adaptation 
  The function then returns the number of cells refined. */

  astats s = adapt_wavelet ({lim}, (double[]){ETAE},
			    MAXLEVEL, MINLEVEL);
  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
  return s.nf;
}


event init (i = 0){
  foreach()
  zb[]=-10;
  conserve_elevation();

  /**
  The initial still water surface is at $z=0$ so that the water depth
  $h$ is... */

  foreach()
    h[] = max(0., - zb[]);
  boundary ({h});

  }

/**
## Outputs

### At each timestep

We output simple summary statistics for *h* and *u.x* on standard
error. */

event logfile (i++) {
  stats s = statsf (h);
  norm n = normf (u.x);
  if (i == 0)
    fprintf (stderr, "t i h.min h.max h.sum u.x.rms u.x.max dt\n");
  fprintf (stderr, "%g %d %g %g %g %g %g %g\n", t, i, s.min, s.max, s.sum, 
	   n.rms, n.max, dt);

}

/**
### Snapshots
*/

event snapshots (i++; i <= 15) {
  /**
  We save snapshot files along the way. */

  char *outfile = NULL;
  outfile = (char *) malloc(sizeof(char) * 16);
  sprintf(outfile, "snapshot-%d.gfs", i);
  output_gfs (file = outfile, t = t);
}

/**
## Adaptivity

And finally we apply our *adapt()* function at every timestep. */

event do_adapt (i++) adapt();
