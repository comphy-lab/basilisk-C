/**
# Initialization of a sphere in a tri-periodic box using MPI with nproc=4

It seems that there is a problem (or is this a feature?) with the
distribution of processes that leads to the simulation being stuck in
the event *init* (more specifically the issue seems to come after
using the adapt_wavelet function), when using 4 procs and a minlevel=1
in adapt_wavelet.

The issue is resolved when using a minlevel>=4 in adapt_wavelet, or
when removing the periodicity in at least one direction, or when using
a number of procs $2^n$, where n is a multiple of 3, or when changing
the position of the object to the opposite side of the domain
(-L0/4. for all dimensions). */

#include "grid/octree.h"
#include "../myembed.h"
#include "../mycentered.h"
#include "view.h"

/**
## Reference solution */

#define d (1.)

/**
We also define the shape of the domain. */

#define sphere(x,y,z) (sq ((x)) + sq ((y)) + sq ((z)) - sq ((d)/2.))
coord p_p;

void p_shape (scalar c, face vector f, coord p)
{
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (sphere ((x - p.x),
		     (y - p.y),
		     (z - p.z)));
  boundary ({phi});
  fractions (phi, c, f);
  fractions_cleanup (c, f,
		     cmin = 1.e-14, smin = 1.e-14);
}

/**
## Setup

We define the mesh adaptation parameters. */

#define lmin (5) // Min mesh refinement level (l=5 is 2pt/d)
#define lmax (7) // Max mesh refinement level (l=7 is 8pt/d)

int main ()
{  
  /**
  The domain is $16\times 16\times 16$ and tri-periodic. */

  L0 = 16.;
  size (L0);
  origin (-(L0)/2., -(L0)/2., -(L0)/2.);
  
  /**
  We initialize the grid. */
  
  N = 1 << (lmin);
  init_grid (N);

  foreach_dimension()
    periodic (left);

  run();
}

/**
## Initial conditions */

event init (i = 0)
{
  /**
  We initialize the embedded boundary. */

  /**
  We first define the sphere's initial position. */

  p_p.x = (L0)/4.;
  p_p.y = (L0)/4.;
  p_p.z = (L0)/4.;

  fprintf (stderr, "%d %g %g %g %g %g %g", i, t, L0, (d), p_p.x, p_p.y, p_p.z), fflush (stderr);
  
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
    
    scalar pid[];
    foreach()
      pid[] = pid();

    clear();
    view (fov = 25, camera = "front",
	  tx = 0., ty = 0.,
	  bg = {1,1,1},
	  width = 800, height = 800);

    draw_vof ("cs", "fs", lw = 5);
    cells (n = {0,0,1}, alpha = (p_p.z));
    squares ("pid", n = {0,0,1}, alpha = (p_p.z));
    save ("mesh-xy.png");

    clear ();
    view (fov = 25, camera = "top",
	  tx = 0., ty = 0.,
	  bg = {1,1,1},
	  width = 800, height = 800);

    draw_vof ("cs", "fs", lw = 5);
    cells (n = {0,1,0}, alpha = (p_p.y));
    squares ("pid", n = {0,1,0}, alpha = (p_p.y));
    save ("mesh-xz.png");

    clear ();
    view (fov = 25, camera = "right",
	  tx = 0., ty = 0.,
	  bg = {1,1,1},
	  width = 800, height = 800);

    draw_vof ("cs", "fs", lw = 5);
    cells (n = {1,0,0}, alpha = (p_p.x));
    squares ("pid", n = {1,0,0}, alpha = (p_p.x));
    save ("mesh-yz.png");
    
  } while ((ss.nf || ss.nc) && ic < 100);
#endif // TREE
  
  p_shape (cs, fs, p_p);
}
