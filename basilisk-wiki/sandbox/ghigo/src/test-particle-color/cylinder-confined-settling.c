/**
# Settling cylinder in an long channel for $r=1.5$

This test case is based on the numerical work of [Yu et al.,
2002](#yu2002) and [Wachs, 2009](#wachs2009). We investigate the
settling of a cylinder of diameter $d$ in a long channel of aspect
ration $W/d = 4$.

This test case is governed by the density ratio $r=\rho_s/\rho$ and
the Reynolds number $Re = \frac{UD}{\nu}$, where $U = \sqrt{\frac{\pi
D}{2}(\frac{\rho_s}{\rho} - 1)g}$ is the "steady" settling velocity.

Due to the added-mass effect for density ratios close to 1, we choose
to simulate the cases where $r=[1.5,3]$, leading to a Reynolds number
$Re=[346.8,693.6]$.

We solve here the 2D Navier-Stokes equations and describe the cylinder
using an [embedded boundary](/src/embed.h). */

/**
## Notes

A minimum level *lmax=12* seems to be necessary to prevent the
particle from crashing into the wall. */

#include "grid/quadtree.h"
#include "../myembed.h"
#include "../mycentered.h"
#include "../myembed-particle-color.h"
#include "../myperfs.h"
#include "view.h"

/**
## Reference solution */

#define d    (1.)
#define h    (4.*(d)) // Width of the channel
#define nu   (0.0080880)
#define grav (10.)
#if DENSITY == 1
#define r    (3.) // Ratio of solid to fluid density
#else
#define r    (1.5) // Ratio of solid to fluid density
#endif // DENSITY
#define uref (sqrt (M_PI*(d)/2.*((r) - 1.)*(grav))) // Characteristic speed
#define tref ((d)/(uref)) // Characteristic time

/**
We also define the shape of the domain. */

#define cylinder(x,y) (sq ((x)) + sq ((y)) - sq ((d)/2.))
#define wall(x,w)     ((x) - (w)) // + over, - under

void p_shape (scalar c, face vector f, coord p)
{
  vertex scalar phi[];
  foreach_vertex() {
    phi[] = HUGE;
    for (double xp = -(L0); xp <= (L0); xp += (L0))
      for (double yp = -(L0); yp <= (L0); yp += (L0))
	phi[] = intersection (phi[],
			      intersection ((cylinder ((x + xp - p.x),
						       (y + yp - p.y))),
					    intersection (
							  -(wall (x,  (h)/2.)),
							   (wall (x, -(h)/2.)))));
  }
  boundary ({phi});
  fractions (phi, c, f);
  fractions_cleanup (c, f,
		     smin = 1.e-14, cmin = 1.e-14);
}

/**
We finally define the particle parameters. */

const double p_r = (r); // Ratio of solid and fluid density
const double p_v = (p_volume_cylinder ((d))); // Particle volume
const coord  p_i = {(p_moment_inertia_cylinder ((d), (r))),
		    (p_moment_inertia_cylinder ((d), (r)))}; // Particle moment of interia
const coord  p_g = {0., -(grav)}; // Gravity, zero

/**
## Setup

We need a field for viscosity so that the embedded boundary metric can
be taken into account. */

face vector muv[];

/**
Since both the walls of the channel and the cylinder are described
using the same embedded volume fraction *cs*, we use the color field
*p_col* to color the cylinder only. */

void p_shape_col (scalar c, coord p)
{
  vertex scalar phi[];
  foreach_vertex() {
    phi[] = HUGE;
    for (double xp = -(L0); xp <= (L0); xp += (L0))
      for (double yp = -(L0); yp <= (L0); yp += (L0))
	phi[] = intersection (phi[],
			      (cylinder ((x + xp - p.x),
					 (y + yp - p.y))));
  }
  boundary ({phi});
  fractions (phi, c);
}

/**
Finally, we define the mesh adaptation parameters. */

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
  As we are computing an equilibrium solution, we remove the Neumann
  pressure boundary condition which is responsible for
  instabilities. */

  for (scalar s in {p, pf}) {
    s.neumann_zero = true;
  }
  
  /**
  We initialize the embedded boundary. */

  /**
  We first define the particle's initial position. */

  p_p.x = -1.;
  p_p.y = (L0) - 10.*(d);
  
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

event logfile (i++; t < 150.*(tref))
{
  coord Fp, Fmu;
  p_shape_col (p_col, (p_p));
  boundary ((scalar *) {p_col});
  embed_color_force (p, u, mu, p_col, &Fp, &Fmu);

  double CD = (Fp.x + Fmu.x)/(0.5*sq ((p_u.y) + SEPS)*(d));
  double CL = (Fp.y + Fmu.y)/(0.5*sq ((p_u.y) + SEPS)*(d));

  fprintf (stderr, "%d %g %g %d %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	   i, t/(tref), dt/(tref),
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resb, mgp.resa,
	   mgu.resb, mgu.resa,
	   p_p.x, p_p.y,
	   p_u.x/(uref), p_u.y/(uref),
	   p_w.x, p_w.y,
	   CD, CL,
	   fabs (p_u.y)*(d)/nu
	   );
  fflush (stderr);

  double cell_wall = fabs (p_p.y - (d)/2.)/((L0)/(1 << (lmax)));
  if (cell_wall <= 1.)
    return 1; // Stop
}

/**
## Results

#### Results for $Re = 346.8$

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
     'log' u ($14 + 2.):($15 - 78.) w l lc rgb "blue" t "Basilisk, l=13"     
~~~

~~~gnuplot Time evolution of the Reynolds number
set key bottom right
set xlabel "t/(d/u)"
set ylabel "Re"
set xrange [0:150]
set yrange [0:600]
plot 259.5 w l lw 2 lc rgb "black" t "Wachs, 2009, r=1.5",		\
     522   w l lw 2 lc rgb "brown" t "Wachs, 2009, r=3",			\
     'log' u 2:22 w l lc rgb "blue" t "Basilisk, l=13"
~~~

~~~gnuplot Time evolution of the drag coefficient
set key top right
set ylabel "C_D"
set yrange [1:5]
plot 1.785 w l lw 2 lc rgb "black" t "Wachs, 2009, r=1.5",		\
     1.764 w l lw 2 lc rgb "brown" t "Wachs, 2009, r=3",			\
     '< cat log | awk -f ../data/Wachs2009/surface.awk' u 1:2 w l lc rgb "blue" t "Basilisk, l=13"
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
