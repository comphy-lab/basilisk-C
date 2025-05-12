/**
# Rising bubble

This case is almost exactly the same as the one proposed by S. Popinet in 2D. All credit to him !

A three-dimensional bubble is released in a box and raises
under the influence of buoyancy. This test case is similar to the case proposed by
[Hysing et al, 2009](/src/references.bib#hysing2009) (see also [the
FEATFLOW
page](http://www.featflow.de/en/benchmarks/cfdbenchmarking/bubble.html)).

We solve the incompressible, variable-density, Navier--Stokes
equations with interfaces and surface tension. We can used standard or "reduced"
gravity. */


#define REDUCED 1
#define CASE2 1
#define ADAPT 0

#include "grid/multigrid3D.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#if REDUCED
# include "reduced.h"
#endif

#ifndef LEVEL
# define LEVEL 8
#endif

/**
The boundary conditions are slip lateral walls (the default) and
no-slip on the right and left walls. */

u.t[right] = dirichlet(0);
u.t[left]  = dirichlet(0);

/**
We make sure there is no flow through the top and bottom boundary,
otherwise the compatibility condition for the Poisson equation can be
violated. */

uf.n[bottom] = 0.;
uf.n[top] = 0.;

int main() {

  /**
  The domain will span $[0:2]\times[0:1]\times[0:1]$ and will be resolved with
  $256\times 128\times 128$ grid points. */

  size (2);
  origin (0., -0.5, -0.5);
  dimensions(nx=4,ny=2,nz=2);
  init_grid (1 << LEVEL);
  
  /**
  Hysing et al. consider two cases (1 and 2), with the densities, dynamic
  viscosities and surface tension of fluid 1 and 2 given below. */

  rho1 = 1000., mu1 = 10.;
#if CASE2
  rho2 = 1., mu2 = 0.1, f.sigma = 1.96;
#else
  rho2 = 100., mu2 = 1., f.sigma = 24.5;
#endif

  /**
  We reduce the tolerance on the Poisson and viscous solvers to
  improve the accuracy. */
  
  TOLERANCE = 1e-4;
#if REDUCED
  G.x = -0.98;
  Z.x = 1.;
#endif
  run();
}

event init (t = 0) {

  /**
  The bubble is centered on (0.5,0,0) and has a radius of 0.25. */

  fraction (f, sq(x - 0.5) + sq(y) + sq(z) - sq(0.25));
}

/**
We add the acceleration of gravity. */

#if !REDUCED
event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] -= 0.98;
}
#endif


/**
We log the position of the center of mass of the bubble, its velocity
and volume as well as convergence statistics for the multigrid
solvers. */

event extract (i++) {
  double xb = 0., vb = 0., sb = 0.;
  foreach(reduction(+:xb) reduction(+:vb) reduction(+:sb)) {
    double dv = (1. - f[])*dv();
    xb += x*dv;
    vb += u.x[]*dv;
    sb += dv;
  }
  double s = interface_area(f);

  char name[80];
  sprintf (name, "velocity_pos.dat");
  static FILE * fp = fopen (name, "w");
  fprintf (fp,"%g %g %g %g %g %g %g\n", 
	  t, sb, -1., xb/sb, vb/sb, s, dt);
  //putchar ('\n');
  //fflush (stdout);
  fflush (fp);
}

event logfile(i++) {
  printf ("i = %d t = %g\n", i,t);
  //printf ("utc = %g", utc);
  fflush(stdout);
}



event interface (t += 1; t <= 3) {

  /**
  We only output the interface in this function. */
  char name[80];
  sprintf (name, "interface-%f-%d.txt", t,pid());
  FILE* fp = fopen (name, "w");
  output_facets (f, fp);
  fclose (fp);
   
  char named[80];
  sprintf (named, "dump-%g", t);
  dump (named);

}


/**
If gfsview is installed on the system, we can also visualise the
simulation as it proceeds. */

#if 0
event gfsview (i += 10) {
  static FILE * fp = popen("gfsview2D rising.gfv", "w");
  scalar vort[];
  vorticity (u, vort);
  output_gfs (fp);
}
#endif

#if ADAPT
event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){5e-4,5e-4,5e-4}, LEVEL);
}
#endif
