/**
# Two cylinder moving at the same speed as the surrounding inviscid flow

In this test case, both the fluid and the cylinders are moving at the
same speed. The presence of the embedded boundary should not create
any disturbance in the flow.

We solve here the Euler equations and add the cylinder using an
[embedded boundary](/src/embed.h). */

#include "../myembed.h"
#include "../mycentered.h"
#include "../myembed-moving-bicolor.h"
#include "view.h"

/**
## Reference solution */

#define d0   (0.753)
#define d1   (0.553)

#define uref (0.912) // Reference velocity, uref
#define tref ((d1)/(uref)) // Reference time, tref=d/u

/**
We also define the shape of the cylinders. */

#define circle(x,y,r) (sq (x) + sq (y) - sq (r))

double p0_phi (double xc, double yc)
{
  return circle (xc, yc, (d0)/2.);
}
double p1_phi (double xc, double yc)
{
  return circle (xc, yc, (d1)/2.);
}

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
  We set the maximum timestep. */

  DT = 1.e-2*(tref);
 
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
  for (scalar s in {u}) {
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

  We define the cylinders' initial position. */

  for (int i = 0; i < (p_n); i++)
    foreach_dimension() 
      p_p[i].x = -10 + 10.*((double) i);

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
  We initialize the cylinders' velocity. */

  for (int i = 0; i < (p_n); i++)
    p_u[i].x = (uref);
  
  /**
  We initialize the velocity to speed-up convergence. */

  foreach()
    u.x[] = (uref);
  boundary ((scalar *) {u});  
}

/**
## Embedded boundaries 

The cylinders' position are advanced to time $t + \Delta t$. */

event advection_term (i++)
{
  for (int i = 0; i < (p_n); i++)
    foreach_dimension()
      p_p[i].x += p_u[i].x*(dt);
}

/**
We also verify here that the velocity and pressure gradient
boundary conditions are correctly computed. */

event check (i++)
{
  foreach() {
    if (cs[] > 0. && cs[] < 1.) {

      // Normal pointing from fluid to solid
      coord b, n;
      embed_geometry (point, &b, &n);
      
      // Velocity
      bool dirichlet;
      double ub;

      // Particle 1
      if (p0_col[] > 0. && p0_col[] < 1.) {
	ub = u.x.boundary[embed] (point, point, u.x, &dirichlet);
	assert (dirichlet);
	assert (ub - p_u[0].x == 0.);
	ub = u.y.boundary[embed] (point, point, u.y, &dirichlet);
	assert (dirichlet);
	assert (ub == 0.);
      }

      // Particle 2
      else if (p1_col[] > 0. && p1_col[] < 1.) {
	ub = u.x.boundary[embed] (point, point, u.x, &dirichlet);
	assert (dirichlet);
	assert (ub - p_u[1].x == 0.);
	ub = u.y.boundary[embed] (point, point, u.y, &dirichlet);
	assert (dirichlet);
	assert (ub == 0.);
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
  
  assert (normf(e).max < 1.e-10);
}

event snapshot (i += 10)
{
  view (fov = 20, camera = "front",
	tx = 0., ty = 0.,
	bg = {1,1,1},
	width = 800, height = 800);

  draw_vof ("cs", "fs", lw = 5);
  squares ("u.x", min = 0);
  save ("mesh.mp4");
}

/**
## Results

![Mesh](bicylinder-steady-moving/mesh.mp4)(loop)

We plot the time evolution of the error. We observe small variations
of the velocity.

~~~gnuplot Time evolution of the average error
reset
set terminal svg font ",16"
set key top left spacing 1.1
set grid ytics
set xtics 0,1,10
set ytics format "%.0e" 1.e-18,1.e-4,1.e4
set xlabel 't/(d/u)'
set ylabel '||error||_{1}'
set xrange [0:2]
set yrange [1.e-18:1.e3]
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
