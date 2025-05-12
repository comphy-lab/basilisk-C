/**
# Test dump and restore functions */ 

#include "grid/octree.h"
#include "../myembed.h"
#include "../mycentered.h"
#include "../myembed-moving.h"
#include "view.h"

/**
## Reference solution */

#define d (1)

/**
We define the shape of the domain. */

void p_shape (scalar c, face vector f, coord pos)
{
  vertex scalar phi[];
  foreach_vertex()
    phi[] = sq (x - pos.x) + sq (y - pos.y) + sq (z - pos.z) - sq ((d)/2.);
  boundary ({phi});
  fractions (phi, c, f);
  fractions_cleanup (c, f,
		     smin = 1.e-14, cmin = 1.e-14);
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
  for (scalar s in {u, p})
    s.third = false;
#else
  for (scalar s in {u, p})
    s.third = true;
#endif // ORDER2

  /**
  We use a slope-limiter to reduce the errors made in small-cells. */
  
#if SLOPELIMITER
  for (scalar s in {u, p}) {
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
  As we are computing an equilibrium solution, we also remove the
  Neumann pressure boundary condition which is responsible for
  instabilities. */

  for (scalar s in {p}) {
    s.neumann_zero = true;
  }
  
  /**
  We initialize the embedded boundary.

  We define the cylinder's initial position. */

  foreach_dimension()
    p_p.x = 4.3;

  #if TREE
  /**
  When using *TREE*, we refine the mesh around the embedded
  boundary. */
  
  astats ss;
  int ic = 0;
  do {
    ic++;
    p_shape (cs, fs, p_p);
    ss = adapt_wavelet ({cs}, (double[]) {1.e-30},
			maxlevel = (lmax), minlevel = (1));
  } while ((ss.nf || ss.nc) && ic < 100);
#endif // TREE
  
  p_shape (cs, fs, p_p);
  
  /**
  We initialize the particle's velocity. */

  foreach_dimension() { 
    p_u.x = -2.7;
    p_w.x = 0.09;
    p_au.x = -1.33;
    p_aw.x = -0.75;
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
  particle pp = {p_p, p_u, p_w, p_au, p_aw};
  struct p_Dump pp_Dump = {"p_dump", &pp};
  p_dump (pp_Dump);

  fprintf (stderr, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	   p_p.x, p_p.y, p_p.z,
	   p_u.x, p_u.y, p_u.z,
	   p_w.x, p_w.y, p_w.z,
	   p_au.x, p_au.y, p_au.z,
	   p_aw.x, p_aw.y, p_aw.z
	   ), fflush (stderr);

  // Restore
  foreach_dimension() { 
    p_p.x = 0.;
    p_u.x = 0.;
    p_w.x = 0.;
    p_au.x = 0.;
    p_aw.x = 0.;
  }

  fprintf (stderr, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	   p_p.x, p_p.y, p_p.z,
	   p_u.x, p_u.y, p_u.z,
	   p_w.x, p_w.y, p_w.z,
	   p_au.x, p_au.y, p_au.z,
	   p_aw.x, p_aw.y, p_aw.z
	   ), fflush (stderr);

  particle pp_restore;
  p_restore ("p_dump", &pp_restore);

  fprintf (stderr, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	   p_p.x, p_p.y, p_p.z,
	   p_u.x, p_u.y, p_u.z,
	   p_w.x, p_w.y, p_w.z,
	   p_au.x, p_au.y, p_au.z,
	   p_aw.x, p_aw.y, p_aw.z
	   ), fflush (stderr);

  foreach_dimension() {
    assert (p_p.x == 4.3);
    assert (p_u.x == -2.7);
    assert (p_w.x == 0.09);
    assert (p_au.x == -1.33);
    assert (p_aw.x == -0.75);
  }
}
