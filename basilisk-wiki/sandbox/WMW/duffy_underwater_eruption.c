/**
# Buoyant plume

~~~gnuplot Maximum height as a function of time
plot 'WG4_minlevel1_maxlevel9_error03_0_0' u 1:2 w l ''
~~~
*/

#include "saint-venant.h"
double default_sea_level=0.;
#include "cgd_read2D.h"
#include "input.h"

/**
We then define a few useful macros and constants. */

#define MAXLEVEL 9
#define MINLEVEL 1
#define ETAE     1e-3 // error on free surface elevation 

int main()
{
  #if QUADTREE
  // 32^2 grid points to start with
  N = 1 << MINLEVEL;
#else // Cartesian
  // 1024^2 grid points
  N = 1 << MAXLEVEL;
#endif
   size (16.);
  origin (-8.,-8.);
  N = 1 << MAXLEVEL;
  run();
}

int adapt() {
#if QUADTREE
  scalar eta[];
  foreach()
    eta[] = h[] > dry ? h[] + zb[] : 0;
  boundary ({eta});

  astats s = adapt_wavelet ({eta}, (double[]){ETAE},
			    MAXLEVEL, MINLEVEL);
  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
  return s.nf;
#else // Cartesian
  return 0;
#endif
}

event init (i = 0)
{

  foreach(){
    zb[]=-1.;
    h[] = max(0., - zb[]);}
  boundary ({h,zb});

/* Open the file. */
  char *fname = "output0.cgd";
  FILE *fidin=fopen(fname,"r");
  deformation_cgd_read (x=0., y=0., fid = fidin, iterate = adapt );
  fprintf(stderr,"deformation\n");
  fclose ( fidin );

}


Gauge gauges[] = {
  {"WG4_minlevel1_maxlevel9_error03_0_0", 0, 0},
  {NULL}
};

event output (t += 0.1; t <= 21)
  output_gauges (gauges, {eta});

/**
## Adaptivity

And finally we apply our *adapt()* function at every timestep. */

event do_adapt (i++) adapt();
