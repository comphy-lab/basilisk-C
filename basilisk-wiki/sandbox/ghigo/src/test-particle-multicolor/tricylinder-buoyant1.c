/**
# Three buoyant moving cylinders advected by an inviscid flow

In this test case, the cylinders are buoyant and therefore should move
with the fluid. The cylinders are initialized with the same speed as
the surrounding fluid and therefore should not create any disturbance
in the flow.

We solve here the Euler equations and add the cylinders using an
[embedded boundary](/src/embed.h). */

#define p_n (3)

#include "../myembed.h"
#include "../mycentered.h"
#include "../myembed-particle-multicolor.h"
#include "view.h"

/**
## Reference solution */

#define d0   (0.753)
#define d1   (0.853)
#define d2   (0.553)

#define uref (0.912) // Reference velocity, uref
#define tref ((d2)/(uref)) // Reference time, tref=d/u

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
double p2_phi (double xc, double yc)
{
  return circle (xc, yc, (d2)/2.);
}

/**
We finally define the gravity acceleration. */

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

u.n[right] = neumann (0);
u.t[right] = neumann (0);
p[right]   = dirichlet (0);

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

  We define the cylinders' shape. */

  pl[0].phi = p0_phi;
  pl[1].phi = p1_phi;
  pl[2].phi = p2_phi;

  /**
  Next we define the cylinders' physical parameters. */

  // Ratio of solid and fluid density
  pl[0].r = 1.;
  pl[1].r = 1.;
  pl[2].r = 1.;
  // Particle volume
  pl[0].v = (p_volume_cylinder ((d0)));
  pl[1].v = (p_volume_cylinder ((d1)));
  pl[2].v = (p_volume_cylinder ((d2)));
  // Particle moment of inertia
  foreach_dimension() {
    pl[0].i.x = (p_moment_inertia_cylinder ((d0), pl[0].r));
    pl[1].i.x = (p_moment_inertia_cylinder ((d1), pl[1].r));
    pl[2].i.x = (p_moment_inertia_cylinder ((d2), pl[2].r));
  }
  
  /**
  We then define the cylinders' initial position. */

  for (int i = 0; i < (p_n); i++)
    foreach_dimension() 
      pl[i].c.x = -10 + 10.*((double) i);

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
  We initialize the cylinders' velocity. */

  for (int i = 0; i < (p_n); i++)
    pl[i].u.x = (uref);
  
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
  
  fprintf (stderr, "%d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	   i, t/(tref), dt/(tref),
	   normf(e).avg, normf(e).max,
	   normf(ep).avg, normf(ep).max,
	   normf(ef).avg, normf(ef).max,
	   pl[0].c.x, pl[0].c.y,
	   pl[0].u.x, pl[0].u.y,
	   pl[0].w.x, pl[0].w.y,
	   pl[1].c.x, pl[1].c.y,
	   pl[1].u.x, pl[1].u.y,
	   pl[1].w.x, pl[1].w.y,
	   pl[2].c.x, pl[2].c.y,
	   pl[2].u.x, pl[2].u.y,
	   pl[2].w.x, pl[2].w.y
	   );
  fflush (stderr);
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

![Mesh](tricylinder-buoyant1/mesh.mp4)(loop)

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
     ''    u 2:($4) w l lw 2 lc rgb "red"   t 'all cells'
~~~

~~~gnuplot Time evolution of the maximum error
set ylabel '||error||_{inf}'
plot 'log' u 2:($7) w l lw 2 lc rgb "black" t 'cut-cells', \
     ''    u 2:($9) w l lw 2 lc rgb "blue"  t 'full cells', \
     ''    u 2:($5) w l lw 2 lc rgb "red"   t 'all cells'
~~~

~~~gnuplot Time evolution of the particles' trajectories
set key top right
set ytics format "%.0f" -40,2,40
set ylabel 'c'
set yrange [-15:15]
unset logscale
plot 'log' u 2:10 w l lw 2 lc rgb "black"     t 'pl[0].c.x', \
     ''    u 2:11 w l lw 2 lc rgb "blue"      t 'pl[0].c.y', \
     ''    u 2:16 w l lw 2 lc rgb "red"       t 'pl[1].c.x', \
     ''    u 2:17 w l lw 2 lc rgb "sea-green" t 'pl[1].c.y', \
     ''    u 2:22 w l lw 2 lc rgb "orange"    t 'pl[2].c.x', \
     ''    u 2:23 w l lw 2 lc rgb "magenta"   t 'pl[2].c.y'
~~~

~~~gnuplot Time evolution of the particles' translation velocities
set ytics format "%.0e" 1.e-18,1.e-4,1.e4
set ylabel 'u'
set yrange [1.e-18:1.e2]
set logscale y
plot 'log' u 2:12              w l lw 2 lc rgb "black"     t 'pl[0].u.x', \
     ''    u 2:(sqrt($13*$13)) w l lw 2 lc rgb "blue"      t 'pl[0].u.y', \
     ''    u 2:18              w l lw 2 lc rgb "red"       t 'pl[1].u.x', \
     ''    u 2:(sqrt($19*$19)) w l lw 2 lc rgb "sea-green" t 'pl[1].u.y', \
     ''    u 2:24              w l lw 2 lc rgb "orange"    t 'pl[2].u.x', \
     ''    u 2:(sqrt($25*$25)) w l lw 2 lc rgb "magenta"   t 'pl[2].u.y'
~~~

~~~gnuplot Time evolution of the particles' rotational velocities
set ytics format "%.0e" 1.e-18,1.e-4,1.e4
set ylabel 'p_w'
set yrange [1.e-18:1.e2]
set logscale y
plot 'log' u 2:(sqrt($14*$14)) w l lw 2 lc rgb "black"     t 'pl[0].w.x', \
     ''    u 2:(sqrt($15*$15)) w l lw 2 lc rgb "blue"      t 'pl[0].w.y', \
     ''    u 2:(sqrt($20*$20)) w l lw 2 lc rgb "red"       t 'pl[1].w.x', \
     ''    u 2:(sqrt($21*$21)) w l lw 2 lc rgb "sea-green" t 'pl[1].w.y', \
     ''    u 2:(sqrt($26*$26)) w l lw 2 lc rgb "orange"    t 'pl[2].w.x', \
     ''    u 2:(sqrt($27*$27)) w l lw 2 lc rgb "magenta"   t 'pl[2].w.y'
~~~
*/
