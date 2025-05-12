/**
# Test dump and restore functions */ 

#define p_n (2)

#include "grid/octree.h"
#include "../myembed.h"
#include "../mycentered.h"
#include "../myembed-moving-multicolor.h"
#include "view.h"

/**
## Reference solution */

#define d1 (1)
#define d2 (2)

/**
We define the shape of the each particle. */

double p0_phi (double xc, double yc, double zc)
{
  return sq (xc) + sq (yc) + sq (zc) - sq ((d1)/2.);
}

double p1_phi (double xc, double yc, double zc)
{
  return sq (xc) + sq (yc) + sq (zc) - sq ((d2)/2.);
}

/**
## Setup

We need a field for viscosity so that the embedded boundary metric can
be taken into account. */

face vector muv[];

/**
We define the mesh adaptation parameters. */

#define lmin (7) // Min mesh refinement level (l=7 is 2pt/D)
#define lmax (9) // Max mesh refinement level (l=9 is 8pt/D)
#define cmax (1.e-2) // Mesh adaptation criterium for velocity

int main ()
{  
  /**
  The domain is $16^3$. */

  L0 = 16.;
  size (L0);
  origin (-L0/2., -L0/2., -L0/2.);

  /**
  We set periodic boundary conditions on all faces. */

  foreach_dimension()
    periodic (left);
  
  /**
  We set the maximum timestep. */

  DT = 1.e-2;
  
  /**
  We set the tolerance of the Poisson solver. */

  TOLERANCE    = 1.e-4;
  TOLERANCE_MU = 1.e-4;
  
  /**
  We initialize the grid. */
  
  N = 1 << (lmin);
  init_grid (N);
  
  run();
}

/**
## Boundary conditions */

/**
## Properties */

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[];
  boundary ((scalar *) {muv});
}

/**
## Initial conditions */

event init (t = 0)
{
  /**
  We set the viscosity field in the event *properties*. */

  mu = muv;

  /**
  We use "third-order" [face flux interpolation](/src/embed.h). */
    
#if ORDER2
  for (scalar s in {u, p, pf})
    s.third = false;
#else
  for (scalar s in {u, p, pf})
    s.third = true;
#endif // ORDER2

  /**
  We use a slope-limiter to reduce the errors made in small-cells. */
  
#if SLOPELIMITER
  for (scalar s in {u, p, pf}) {
    s.gradient = minmod2;
  }
#endif // SLOPELIMITER
    
#if TREE
  /**
  When using *TREE* and in the presence of embedded boundaries, we
  should also define the gradient of *u* at the cell center of
  cut-cells. */
#endif // TREE
  
  /**
  We initialize the embedded boundary.

  We define the cylinder's shape. */

  pl[0].phi = p0_phi;
  pl[1].phi = p1_phi;

  /**
  Next we define the cylinder's initial position. */

  for (int i = 0; i < (p_n); i++)
    foreach_dimension() 
      pl[i].c.x = 4.3 + 2.*i;

  #if TREE
  /**
  When using *TREE*, we refine the mesh around the embedded
  boundary. */
  
  astats ss;
  int ic = 0;
  do {
    ic++;
    p_shape (cs, fs, pl);
    ss = adapt_wavelet ({cs}, (double[]) {1.e-30},
			maxlevel = (lmax), minlevel = (1));
  } while ((ss.nf || ss.nc) && ic < 100);
#endif // TREE
  
  p_shape (cs, fs, pl);

  /**
  We initialize the particle's velocity. */
  
  for (int i = 0; i < (p_n); i++)
    foreach_dimension() { 
      pl[i].u.x  = -2.7 + 2.*i;
      pl[i].w.x  = 0.09 + 2.*i;
      pl[i].au.x = -1.33 + 2.*i;
      pl[i].aw.x = -0.75 + 2.*i;
  }
}

/**
## Embedded boundaries */

/**
## Adaptive mesh refinement */

#if TREE
event adapt (i++)
{
  adapt_wavelet ({cs,u}, (double[]) {1.e-2,(cmax),(cmax)},
  		 maxlevel = (lmax), minlevel = (1));

  /**
  We do not need here to reset the embedded fractions to avoid
  interpolation errors on the geometry as the is already done when
  moving the embedded boundaries. It might be necessary to do this
  however if surface forces are computed around the embedded
  boundaries. */
}
#endif // TREE

/**
## Outputs */

event logfile (i = 1)
{
  // Dump fluid
  dump ();

  // Dump particle
  struct p_Dump pp_Dump = {"p_dump", pl};
  p_dump (pp_Dump);

  fprintf (stderr, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	   pl[0].c.x, pl[0].c.y, pl[0].c.z,
	   pl[1].c.x, pl[1].c.y, pl[1].c.z,
	   pl[0].u.x, pl[0].u.y, pl[0].u.z,
	   pl[1].u.x, pl[1].u.y, pl[1].u.z,	   
	   pl[0].w.x, pl[0].w.y, pl[0].w.z,
	   pl[1].w.x, pl[1].w.y, pl[1].w.z,	   
	   pl[0].au.x, pl[0].au.y, pl[0].au.z,
	   pl[1].au.x, pl[1].au.y, pl[1].au.z,	   
	   pl[0].aw.x, pl[0].aw.y, pl[0].aw.z,
	   pl[1].aw.x, pl[1].aw.y, pl[1].aw.z	   
	   ), fflush (stderr);

  for (int i = 0; i < (p_n); i++)
    foreach_dimension() {
      pl[i].c.x  = 0.;
      pl[i].u.x  = 0.;
      pl[i].w.x  = 0.;
      pl[i].au.x = 0.;
      pl[i].aw.x = 0.;
    }

  fprintf (stderr, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	   pl[0].c.x, pl[0].c.y, pl[0].c.z,
	   pl[1].c.x, pl[1].c.y, pl[1].c.z,
	   pl[0].u.x, pl[0].u.y, pl[0].u.z,
	   pl[1].u.x, pl[1].u.y, pl[1].u.z,	   
	   pl[0].w.x, pl[0].w.y, pl[0].w.z,
	   pl[1].w.x, pl[1].w.y, pl[1].w.z,	   
	   pl[0].au.x, pl[0].au.y, pl[0].au.z,
	   pl[1].au.x, pl[1].au.y, pl[1].au.z,	   
	   pl[0].aw.x, pl[0].aw.y, pl[0].aw.z,
	   pl[1].aw.x, pl[1].aw.y, pl[1].aw.z	   
	   ), fflush (stderr);

  p_restore ("p_dump", pl);

  fprintf (stderr, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	   pl[0].c.x, pl[0].c.y, pl[0].c.z,
	   pl[1].c.x, pl[1].c.y, pl[1].c.z,
	   pl[0].u.x, pl[0].u.y, pl[0].u.z,
	   pl[1].u.x, pl[1].u.y, pl[1].u.z,	   
	   pl[0].w.x, pl[0].w.y, pl[0].w.z,
	   pl[1].w.x, pl[1].w.y, pl[1].w.z,	   
	   pl[0].au.x, pl[0].au.y, pl[0].au.z,
	   pl[1].au.x, pl[1].au.y, pl[1].au.z,	   
	   pl[0].aw.x, pl[0].aw.y, pl[0].aw.z,
	   pl[1].aw.x, pl[1].aw.y, pl[1].aw.z	   
	   ), fflush (stderr);

  for (int i = 0; i < (p_n); i++)
    foreach_dimension() {
      assert (pl[i].c.x  == 4.3 + 2.*i);
      assert (pl[i].u.x  == -2.7 + 2.*i);
      assert (pl[i].w.x  == 0.09 + 2.*i);
      assert (pl[i].au.x == -1.33 + 2.*i);
      assert (pl[i].aw.x == -0.75 + 2.*i);
    }
}
