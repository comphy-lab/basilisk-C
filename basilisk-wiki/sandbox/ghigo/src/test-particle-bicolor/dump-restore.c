/**
# Test dump and restore functions */ 

#include "grid/octree.h"
#include "../myembed.h"
#include "../mycentered.h"
#include "../myembed-moving-bicolor.h"
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

  We define the cylinder's initial position. */

  for (int i = 0; i < (p_n); i++)
    foreach_dimension() 
      p_p[i].x = 4.3 + 2.*i;

  #if TREE
  /**
  When using *TREE*, we refine the mesh around the embedded
  boundary. */
  
  astats ss;
  int ic = 0;
  do {
    ic++;
    p_shape (cs, fs);
    ss = adapt_wavelet ({cs}, (double[]) {1.e-30},
			maxlevel = (lmax), minlevel = (1));
  } while ((ss.nf || ss.nc) && ic < 100);
#endif // TREE
  
  p_shape (cs, fs);
  
  /**
  We initialize the particle's velocity. */

  for (int i = 0; i < (p_n); i++)
    foreach_dimension() { 
      p_u[i].x = -2.7 + 2.*i;
      p_w[i].x = 0.09 + 2.*i;
      p_au[i].x = -1.33 + 2.*i;
      p_aw[i].x = -0.75 + 2.*i;
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
  particle pp[p_n];
  for (int i = 0; i < (p_n); i++) {
    pp[i].c = p_p[i];
    pp[i].u = p_u[i];
    pp[i].w = p_w[i];
    pp[i].au = p_au[i];
    pp[i].aw = p_aw[i];
  }

  struct p_Dump pp_Dump = {"p_dump", pp};
  p_dump (pp_Dump);

  fprintf (stderr, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	   p_p[0].x, p_p[0].y, p_p[0].z,
	   p_p[1].x, p_p[1].y, p_p[1].z,
	   p_u[0].x, p_u[0].y, p_u[0].z,
	   p_u[1].x, p_u[1].y, p_u[1].z,	   
	   p_w[0].x, p_w[0].y, p_w[0].z,
	   p_w[1].x, p_w[1].y, p_w[1].z,	   
	   p_au[0].x, p_au[0].y, p_au[0].z,
	   p_au[1].x, p_au[1].y, p_au[1].z,	   
	   p_aw[0].x, p_aw[0].y, p_aw[0].z,
	   p_aw[1].x, p_aw[1].y, p_aw[1].z	   
	   ), fflush (stderr);

  for (int i = 0; i < (p_n); i++)
    foreach_dimension() {
      p_p[i].x  = 0.;
      p_u[i].x  = 0.;
      p_w[i].x  = 0.;
      p_au[i].x = 0.;
      p_aw[i].x = 0.;
    }

  fprintf (stderr, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	   p_p[0].x, p_p[0].y, p_p[0].z,
	   p_p[1].x, p_p[1].y, p_p[1].z,
	   p_u[0].x, p_u[0].y, p_u[0].z,
	   p_u[1].x, p_u[1].y, p_u[1].z,	   
	   p_w[0].x, p_w[0].y, p_w[0].z,
	   p_w[1].x, p_w[1].y, p_w[1].z,	   
	   p_au[0].x, p_au[0].y, p_au[0].z,
	   p_au[1].x, p_au[1].y, p_au[1].z,	   
	   p_aw[0].x, p_aw[0].y, p_aw[0].z,
	   p_aw[1].x, p_aw[1].y, p_aw[1].z	   
	   ), fflush (stderr);

  particle pp_restore;
  p_restore ("p_dump", &pp_restore);

  fprintf (stderr, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	   p_p[0].x, p_p[0].y, p_p[0].z,
	   p_p[1].x, p_p[1].y, p_p[1].z,
	   p_u[0].x, p_u[0].y, p_u[0].z,
	   p_u[1].x, p_u[1].y, p_u[1].z,	   
	   p_w[0].x, p_w[0].y, p_w[0].z,
	   p_w[1].x, p_w[1].y, p_w[1].z,	   
	   p_au[0].x, p_au[0].y, p_au[0].z,
	   p_au[1].x, p_au[1].y, p_au[1].z,	   
	   p_aw[0].x, p_aw[0].y, p_aw[0].z,
	   p_aw[1].x, p_aw[1].y, p_aw[1].z	   
	   ), fflush (stderr);

  for (int i = 0; i < (p_n); i++)
    foreach_dimension() {
      assert (p_p[i].x  == 4.3 + 2.*i);
      assert (p_u[i].x  == -2.7 + 2.*i);
      assert (p_w[i].x  == 0.09 + 2.*i);
      assert (p_au[i].x == -1.33 + 2.*i);
      assert (p_aw[i].x == -0.75 + 2.*i);
    }
}
