/**
# .cgd read test file.

This works in serial currently but not in parallel
parallel compiling that doesn't work:
qcc -g -O2 -fopenmp -Wall cgd_read_parallel_test.c -o cgd_read_parallel_test.exe -lm
serial compiling that does
qcc -g -O2 -Wall cgd_read_parallel_test.c -o cgd_read_parallel_test.exe -lm

*/
#include "saint-venant.h"
#include "cgd_read.h"
#include "input.h"

/**
We then define a few useful macros and constants. */

#define MAXLEVEL 8
#define MINLEVEL 4
#define ETAE     1e-2 // error on free surface elevation (1 cm)

int main()
{
 size (1.);
  origin (0.,0.);
#if QUADTREE
  // 32^2 grid points to start with
  init_grid (1 << MINLEVEL);
#else // Cartesian
  // 1024^2 grid points
  init_grid (1 << MAXLEVEL);
#endif

  run();
}


int adapt() {
#if QUADTREE
 
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
  char *fname = "/data/basilisk_test/test.cgd";
  FILE *fidin=fopen(fname,"r");
  fprintf(stderr,"deformation\n");
  deformation_cgd_read (x=0., y=0., fid = fidin, iterate = adapt );
  fclose ( fidin );

}

event logfile (i=0) {
  stats s = statsf (h);
  norm n = normf (u.x);
  if (i == 0)
    fprintf (stderr, "t i h.min h.max h.sum u.x.rms u.x.max dt\n");
  fprintf (stderr, "%g %d %g %g %g %g %g %g\n", t, i, s.min, s.max, s.sum, 
	   n.rms, n.max, dt);
}


/**
## Adaptivity

And finally we apply our *adapt()* function at every timestep. */

event do_adapt (i++) adapt();
