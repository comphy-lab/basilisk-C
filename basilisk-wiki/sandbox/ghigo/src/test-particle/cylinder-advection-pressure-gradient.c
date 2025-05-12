/**
# Heavy cylinder advected by a pressure gradient for $Re=20$

In this test case, the cylinder is twice as heavy as the fluid and is
advected in the xy-direction by a pressure gradient.

We solve here the 2D Navier-Stokes equations and describe the cylinder
using an [embedded boundary](/src/embed.h). */ 

#include "../myembed.h"
#include "../mycentered.h"
#include "../myembed-particle.h"
#include "view.h"

/**
## Reference solution */

#define d    (0.753)
#define Re   (20.)
#define uref (0.912) // Reference velocity, uref
#define tref ((d)/(uref)) // Reference time, tref=d/u

/**
We also define the shape of the domain. */

#define cylinder(x,y) (sq ((x)) + sq ((y)) - sq ((d)/2.))

void p_shape (scalar c, face vector f, coord p)
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
  fractions (phi, c, f);
  fractions_cleanup (c, f,
		     smin = 1.e-14, cmin = 1.e-14);
}

/**
We finally define the particle parameters. */

const double p_r = (2.); // Ratio of solid and fluid density
const double p_v = (p_volume_cylinder ((d))); // Particle volume
const coord  p_i = {(p_moment_inertia_cylinder ((d), 2.)),
		    (p_moment_inertia_cylinder ((d), 2.))}; // Particle moment of interia
const coord  p_g = {0., 0.}; // Gravity, zero

/**
## Setup

We need a field for viscosity so that the embedded boundary metric can
be taken into account. */

face vector muv[];

/**
We define the mesh adaptation parameters. */

#define lmin (5) // Min mesh refinement level (l=5 is 3pt/d)
#define lmax (8) // Max mesh refinement level (l=8 is 24pt/d)
#define cmax (1.e-2*(uref)) // Absolute refinement criteria for the velocity field

int main ()
{  
  /**
  The domain is $8\times 8$ and periodic. */

  L0 = 8.;
  size (L0);
  origin (-L0/2., -L0/2.);

  foreach_dimension()
    periodic (left);

  /**
  We set the maximum timestep. */

  DT = 1.e-2*(tref);
 
  /**
  We set the tolerance of the Poisson solver. It seems that reducing
  the tolerance for the viscous solver avoids creating small spurious
  jumps in the rotational velocity. */

  TOLERANCE    = 1.e-6;
  TOLERANCE_MU = 1.e-6*(uref);
  
  /**
  We initialize the grid. */
  
  N = 1 << (lmax);
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
    muv.x[] = (uref)*(d)/(Re)*fm.x[];
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
  We set the acceleration vector. */

  const face vector av[] = {(uref)/sqrt (2.)/(L0),
			    (uref)/sqrt (2.)/(L0)};
  a = av;

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
  We remove the Neumann pressure boundary condition which is
  responsible for instabilities. */

  for (scalar s in {p, pf}) {
    s.neumann_zero = true;
  }
  
  /**
  We initialize the embedded boundary. */
  
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

event logfile (i++; t < 20.*(tref))
{
  double nu = sqrt (sq (p_u.x) + sq (p_u.y));
  nu /= ((uref)*(d)/(Re))/(d);
  fprintf (stderr, "%d %g %g %d %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g %g\n",
	   i, t/(tref), dt/(tref),
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resb, mgp.resa,
	   mgu.resb, mgu.resa,
	   p_p.x, p_p.y,
	   p_u.x/(uref), p_u.y/(uref),
	   p_w.x, p_w.y,
	   nu
	   );
  fflush (stderr);
}

/**
## Results

~~~gnuplot Time evolution of the particle's position
reset
set terminal svg font ",16"
set key bottom right spacing 1.1
set grid ytics
set xtics 0,2,20
set xlabel 't/(d/u)'
set ylabel '{x, y}'
set xrange [0:20]
plot 'log' u 2:14 w l lw 2 lc rgb "black" t 'x', \
     ''    u 2:15 w l lw 2 lc rgb "blue"  t 'y'
~~~

~~~gnuplot Time evolution of the particle's velocity
set ylabel '{u_{p,x}, u_{p,y}}'
plot 'log' u 2:16 w l lw 2 lc rgb "black" t 'u_{p,x}', \
     ''    u 2:17 w l lw 2 lc rgb "blue"  t 'u_{p,y}'
~~~

~~~gnuplot Time evolution of the particle's rotation rate
set ylabel '{w_x, w_y} = w_z'
plot 'log' u 2:18 w l lw 2 lc rgb "black" t 'w_x', \
     ''    u 2:19 w l lw 2 lc rgb "blue"  t 'w_y'
~~~

~~~gnuplot Time evolution of the norm of the particle's velocity
set ylabel '||u_p||_{2}'
plot 'log' u 2:20 w l lw 2 lc rgb "black" notitle
~~~
*/
