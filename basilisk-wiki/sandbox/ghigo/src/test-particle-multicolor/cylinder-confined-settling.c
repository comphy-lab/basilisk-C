/**
# Settling cylinder in an long channel for $r=1.5$

This test case is based on the numerical work of [Yu et al.,
2002](#yu2002) and [Wachs, 2009](#wachs2009). We investigate the
settling of a cylinder of diameter $d$ in a long channel of aspect
ration $w/d = 4$.

This test case is governed by the density ratio $r=\rho_s/\rho$ and
the Reynolds number $Re = \frac{u_{\mathrm{ref}}d}{\nu}$, where
$u_{\mathrm{ref}} = \sqrt{\frac{\pi d}{2}(\frac{\rho_s}{\rho} - 1)g}$
is the "steady" settling velocity.

Due to the added-mass effect for density ratios close to 1, we choose
to simulate the cases where $r=[1.5,3]$, leading to a Reynolds number
$Re=[346.8,693.6]$.

We solve here the 2D Navier-Stokes equations and describe the cylinder
using an [embedded boundary](/src/embed.h). */

/**
## Notes

A minimum level *lmax=12* seems to be necessary to prevent the
particle from crashing into the wall. */

#define p_n (2)

#include "grid/quadtree.h"
#include "../myembed.h"
#include "../mycentered.h"
#include "../myembed-particle-multicolor.h"
#include "../myperfs.h"
#include "view.h"

/**
## Reference solution */

#define d    (1.)
#define h    (4.*(d)) // Width of the channel
#define nu   (0.0080880)
#define grav (10.)
#if DENSITY == 1
#define rr   (3.) // Ratio of solid to fluid density
#else
#define rr   (1.5) // Ratio of solid to fluid density
#endif // DENSITY
#define uref (sqrt (M_PI*(d)/2.*((rr) - 1.)*(grav))) // Characteristic speed
#define tref ((d)/(uref)) // Characteristic time

/**
We also define the shape of the domain. The first particle is in fact
the channel, while the second particle is the cylinder. */

#define circle(x,y) (sq ((x)) + sq ((y)) - sq ((d)/2.))
#define wall(x,w)   ((x) - (w)) // + over, - under

// Channel
double p0_phi (double xc, double yc)
{
  return intersection (-(wall (xc,  (h)/2.)), (wall (xc, -(h)/2.)));
}
// Cylinder
double p1_phi (double xc, double yc)
{
  return circle (xc, yc);
}

/**
We finally define the gravity acceleration. */

const coord  p_g = {0., -(grav)}; // Gravity

/**
## Setup

We need a field for viscosity so that the embedded boundary metric can
be taken into account. */

face vector muv[];

/**
We define the mesh adaptation parameters. */

#define lmin (8)  // Min mesh refinement level (l=8 is 2pt/d)
#define lmax (13) // Max mesh refinement level (l=13 is 64pt/d)
#define cmax (5.e-3*(uref)) // Absolute refinement criteria for the velocity field

int main ()
{
  /**
  The domain is $128\times 128$. */

  L0 = 128.;
  size (L0);
  origin (-L0/2., 0.);

  /**
  We set the maximum timestep. */

  DT = 1.e-2*(tref);
 
  /**
  We set the tolerance of the Poisson solver. */

  TOLERANCE    = 1.e-4;
  TOLERANCE_MU = 1.e-4*(uref);
  
  /**
  We initialize the grid. */

  N = 1 << (lmin);
  init_grid (N);
    
  run();
}

/**
## Boundary conditions

We use no-slip boundary conditions. */

u.n[bottom] = dirichlet (0);
u.t[bottom] = dirichlet (0);

u.n[top] = dirichlet (0);
u.t[top] = dirichlet (0);

/**
We give boundary conditions for the face velocity to "potentially"
improve the convergence of the multigrid Poisson solver. */

uf.n[bottom] = 0;
uf.n[top]    = 0;

/**
## Properties */

event properties (i++)
{
  foreach_face()
    muv.x[] = (nu)*fm.x[];
  boundary ((scalar *) {muv});
}

/**
## Initial conditions */

event init (i = 0)
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

  We define the shape of the two particles. */

  pl[0].phi = p0_phi;
  pl[1].phi = p1_phi;

  /**
  Next we define the particles' physical parameters. The channel is
  defined using the default parameters. */

  // Ratio of solid and fluid density
  pl[1].r = (rr);
  // Particle volume
  pl[1].v = (p_volume_cylinder ((d)));
  // Particle moment of inertia
  foreach_dimension()
    pl[1].i.x = (p_moment_inertia_cylinder ((d), pl[1].r));

  /**
  We then define the cylinder's initial position. */

  pl[1].c.x = -1.;
  pl[1].c.y = (L0) - 10.*(d);
  
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
}

/**
## Embedded boundaries 

We verify here that the velocity and pressure gradient boundary
conditions for the channel particle are correctly set. */

event check (i++)
{
  foreach() {
    if (cs[] > 0. && cs[] < 1.) {

      // Velocity
      bool dirichlet;
      double ub;

      scalar col = pl[0].col;
      if (col[] > 0. && col[] < 1.) {
	ub = u.x.boundary[embed] (point, point, u.x, &dirichlet);
	assert (dirichlet);
	assert (ub == 0. && pl[0].u.x == 0.);
	ub = u.y.boundary[embed] (point, point, u.y, &dirichlet);
	assert (dirichlet);
	assert (ub == 0. && pl[0].u.y == 0.);
      }

      // Pressure
      bool neumann;
      double pb;
      
      pb = p.boundary[embed] (point, point, p, &neumann);
      assert (!neumann);
      assert (pb == 0.);
      
      // Pressure gradient
      double gb;
      
      gb = g.x.boundary[embed] (point, point, g.x, &dirichlet);
      assert (dirichlet);
      assert (gb == 0.);
      gb = g.y.boundary[embed] (point, point, g.y, &dirichlet);
      assert (dirichlet);
      assert (gb == 0.);
    }
  }
}

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

event logfile (i++; t < 150.*(tref))
{
  p_shape_col (pl);

  coord Fp, Fmu;
  embed_color_force (p, u, mu, pl[1].col, &Fp, &Fmu);

  double CD = (Fp.x + Fmu.x)/(0.5*sq ((pl[1].u.y) + SEPS)*(d));
  double CL = (Fp.y + Fmu.y)/(0.5*sq ((pl[1].u.y) + SEPS)*(d));

  fprintf (stderr, "%d %g %g %d %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	   i, t/(tref), dt/(tref),
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resb, mgp.resa,
	   mgu.resb, mgu.resa,
	   pl[1].c.x, pl[1].c.y,
	   pl[1].u.x/(uref), pl[1].u.y/(uref),
	   pl[1].w.x, pl[1].w.y,
	   CD, CL,
	   fabs (pl[1].u.y)*(d)/nu
	   );
  fflush (stderr);

  double cell_wall = fabs (pl[1].c.y - (d)/2.)/((L0)/(1 << (lmax)));
  if (cell_wall <= 1.)
    return 1; // Stop
}

event animation (i += 10)
{
  scalar omega[];
  vorticity (u, omega);
  boundary ({omega});
  
  view (fov = 2, camera = "front",
	tx = 0., ty = -(pl[1].c.y + 3*(d))/(L0),
	bg = {1,1,1},
	width = 800, height = 800);

  draw_vof ("cs", "fs", filled = -1);
  squares ("omega", map = cool_warm);
  save ("vorticity.mp4");
}

/**
## Results

![Vorticity](cylinder-confined-settling/vorticity.mp4)(loop)

~~~gnuplot Time evolution of the particle's trajectory
reset
set terminal svg font ",16"
set key top right spacing 1.1
set xlabel 'x/d'
set ylabel 'y/d'
set xrange [0.9:3]
set yrange [-60:60]
plot "../data/Wachs2009/Wachs2009-fig5a-r-1p5.csv" u 1:2 w l lw 1 lc rgb "black" t "fig. 5a, Wachs, 2009, r=1.5", \
     "../data/Wachs2009/Wachs2009-fig5a-r-3.csv"   u 1:2 w l lw 1 lc rgb "brown" t "fig. 5a, Wachs, 2009, r=3", \
     "../../test-particle-color/cylinder-confined-settling/log" u ($14 + 2.):($15 - 78.) w l lc rgb "blue" t "Basilisk, l=13", \
     'log' u ($14 + 2.):($15 - 78.) w l lc rgb "red" t "Basilisk, l=13"
~~~

~~~gnuplot Time evolution of the Reynolds number
set key bottom right
set xlabel "t/(d/u)"
set ylabel "Re"
set xrange [0:150]
set yrange [0:600]
plot 259.5 w l lw 2 lc rgb "black" t "Wachs, 2009, r=1.5",		\
     522   w l lw 2 lc rgb "brown" t "Wachs, 2009, r=3", \
     '../../test-particle-color/cylinder-confined-settling/log' u 2:22 w l lc rgb "blue" t "Basilisk, l=13", \
     'log' u 2:22 w l lc rgb "red" t "Basilisk, l=13"
~~~

~~~gnuplot Time evolution of the drag coefficient
set key top right
set ylabel "C_D"
set yrange [1:5]
plot 1.785 w l lw 2 lc rgb "black" t "Wachs, 2009, r=1.5",		\
     1.764 w l lw 2 lc rgb "brown" t "Wachs, 2009, r=3",			\
     '< cat ../../test-particle-color/cylinder-confined-settling/log | awk -f ../data/Wachs2009/surface.awk' u 1:2 w l lc rgb "blue" t "Basilisk, l=13", \
     '< cat log | awk -f ../data/Wachs2009/surface.awk' u 1:2 w l lc rgb "red" t "Basilisk, l=13"
~~~

## References

~~~bib
@article{yu2002,
  title={Viscoelastic mobility problem of a system of particles},
  author={Yu, Z. and P.-T., N. and F., Y. and T., R.},
  journal={Journal of Non-Newtonian Fluid Mechanics},
  volume={104},
  pages={87--124},
  year={2002}
}

@article{wachs2009,
  title={A DEM-DLM/FD method for direct numerical simulation of particulate flows: Sedimentation of polygonal isometric particles in a Newtonian fluid with collisions},
  author={Wachs, A.},
  journal={Computers & Fluids},
  volume={38},
  pages={1608--1628},
  year={2009}
}
*/
