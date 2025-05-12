/**
# Buoyant cylinder moving at the same speed as the surrounding inviscid flow

In this test case, the cylinder is buoyant and therefore should move
with the fluid. The cylinder is initialized with the same speed as the
surrounding fluid and therefore should not create any disturbance in
the flow.

We solve here the Euler equations and add the cylinder using an
[embedded boundary](/src/embed.h). */

#include "../myembed.h"
#include "../mycentered.h"
#include "../myembed-particle.h"
#include "view.h"

/**
## Reference solution */

#define d    (0.753)
#define uref (0.912) // Reference velocity, uref
#define tref ((d)/(uref)) // Reference time, tref=d/u

/**
We also define the shape of the domain. */

#define cylinder(x,y) (sq ((x)) + sq ((y)) - sq ((d)/2.))

void p_shape (scalar c, face vector f, coord p)
{
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (cylinder ((x - p.x), (y - p.y)));
  boundary ({phi});
  fractions (phi, c, f);
  fractions_cleanup (c, f,
		     smin = 1.e-14, cmin = 1.e-14);
}

/**
We finally define the particle parameters. */

const double p_r = (1.); // Ratio of solid and fluid density
const double p_v = (p_volume_cylinder ((d))); // Particle volume
const coord  p_i = {(p_moment_inertia_cylinder ((d), 1.)),
		    (p_moment_inertia_cylinder ((d), 1.))}; // Particle moment of interia
const coord  p_g = {751., 83.6}; // Gravity, random

/**
## Setup

We define the mesh adaptation parameters. */

#define lmin (7) // Min mesh refinement level (l=7 is 3pt/d)
#define lmax (10) // Max mesh refinement level (l=10 is 24pt/d)
#define cmax (1.e-2*(uref)) // Absolute refinement criteria for the velocity field

int main ()
{  
  /**
  The domain is $32\times 32$. */

  L0 = 32.;
  size (L0);
  origin (-L0/2., -L0/2.);

  /**
  We set the maximum timestep. Since we are computing an equilibrium
  solution, we reduce the time step to avoid temporal instabilities
  due to the explicit first-order coupling. */

  DT = 1.e-3*(tref);
 
  /**
  We set the tolerance of the Poisson solver. */

  TOLERANCE    = 1.e-4;
  TOLERANCE_MU = 1.e-4*(uref);
  
  /**
  We initialize the grid. */
  
  N = 1 << (lmax);
  init_grid (N);
  
  run();
}

/**
## Boundary conditions

We use inlet boundary conditions. */

u.n[left] = dirichlet ((uref));
u.t[left] = dirichlet (0);
p[left]   = neumann (0);
pf[left]  = neumann (0);

u.n[right] = neumann (0);
u.t[right] = neumann (0);
p[right]   = dirichlet (0);
pf[right]  = dirichlet (0);

/**
We give boundary conditions for the face velocity to "potentially"
improve the convergence of the multigrid Poisson solver. */

uf.n[left]   = (uref);
uf.n[bottom] = 0;
uf.n[top]    = 0;

/**
## Properties */

/**
## Initial conditions */

event init (i = 0)
{
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
			maxlevel = (lmax), minlevel = (lmin));
  } while ((ss.nf || ss.nc) && ic < 100);
#endif // TREE
  
  p_shape (cs, fs, p_p);

  /**
  We initialize the particle's velocity. */

  p_u.x = (uref);

  /**
  We initialize the velocity to speed-up convergence. */

  foreach()
    u.x[] = (uref);
  boundary ((scalar *) {u});  
}

/**
## Embedded boundaries */

/**
## Adaptive mesh refinement */

#if TREE
event adapt (i++)
{
  adapt_wavelet ({cs,u}, (double[]) {1.e-2,(cmax),(cmax)},
  		 maxlevel = (lmax), minlevel = (lmin));

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

event logfile (i++; t < 2.*(tref))
{
  scalar e[], ef[], ep[];
  foreach() {
    if (cs[] <= 0.)
      e[] = ef[] = ep[] = nodata;
    else {
      e[] = sqrt (sq (u.x[] - (uref)) + sq (u.y[]));
      ep[] = cs[] < 1. ? e[] : nodata;
      ef[] = cs[] >= 1. ? e[] : nodata;
    }
  }
  boundary ((scalar *) {e, ef, ep});
  
  fprintf (stderr, "%d %g %g %g %g %g %g %g %g\n",
	   i, t/(tref), dt/(tref),
	   normf(e).avg, normf(e).max,
	   normf(ep).avg, normf(ep).max,
	   normf(ef).avg, normf(ef).max
	   );
  fflush (stderr);

  /**
  Criteria on maximum value of error. */
  
  assert (normf(e).max < 1.e-8);
}

/**
## Results

We plot the time evolution of the error. We observe small variations
of the velocity, contrary to the equivalent fixed cylinder test case.

~~~gnuplot Time evolution of the average error
reset
set terminal svg font ",16"
set key top right spacing 1.1
set grid ytics
set xtics 0,1,10
set ytics format "%.0e" 1.e-18,1.e-2,1.e-0
set xlabel 't/(d/u)'
set ylabel '||error||_{1}'
set xrange [0:2]
set yrange [1.e-18:1.e-6]
set logscale y
plot 'log' u 2:($6) w l lw 2 lc rgb "black" t 'cut-cells', \
     ''    u 2:($8) w l lw 2 lc rgb "blue"  t 'full cells', \
     ''    u 2:($4) w l lw 2 lc rgb "red"   t 'all cells
~~~

~~~gnuplot Time evolution of the maximum error
set ylabel '||error||_{inf}'
plot 'log' u 2:($7) w l lw 2 lc rgb "black" t 'cut-cells', \
     ''    u 2:($9) w l lw 2 lc rgb "blue"  t 'full cells', \
     ''    u 2:($5) w l lw 2 lc rgb "red"   t 'all cells
~~~
*/
